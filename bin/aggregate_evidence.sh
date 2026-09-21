#!/usr/bin/env bash
# aggregate_evidence.sh — combine geNomad, diamond (U-RVDB-prot, UniRef90),
# protein header taxmaps, bowtie2 coverage, and viroid blastn into one
# evidence table.
#
# Layout: one data row per contig, with one row per U-RVDB-prot hit when hits
# exist and one row per viroid blastn HSP when hits exist. A contig with both
# (rare: the legs target disjoint size classes) gets rvdb rows and viroid
# rows; each row carries NA for the other leg. geNomad, UniRef90, and bowtie2
# columns repeat across rows of the same contig.
#
#   usage: aggregate_evidence.sh -s <sample> -o <outdir> -d <dbdir> [--rvdb-taxmap <file>]
#
# Inputs (default layout):
#   $out/05_genomad/<s>/<s>_virus_summary.tsv
#   $out/06_bowtie2/<s>/<s>.coverage.tsv
#   $out/07_diamond/<s>/<s>.rvdb.tsv
#   $out/07_diamond/<s>/<s>.uniref90.tsv
#   $out/09_viroid/<s>/<s>.viroid.tsv          (optional; 16-field blastn outfmt)
#   $db/blastx/uniref90.taxmap.tsv
#   $db/blastx/U-RVDBv32.0-prot.taxmap.tsv
#   $db/blastn/*_seqids.txt                    (optional; "accession<TAB-or-space>description")
#
# One-time sorted caches built beside the DBs on first run (idempotent):
#   $db/blastx/U-RVDBv32.0-prot.taxmap.srt
#   $db/blastx/uniref90.taxmap.srt
#   $db/blastn/viroid_seqids.srt
#
# Output: $out/evidence/<s>.evidence.tsv  (73-column schema, see header)
set -euo pipefail
export LC_ALL=C
export TMPDIR="$PWD"

usage() { cat <<'EOF'
Usage: aggregate_evidence.sh -s <sample> -o <outdir> -d <dbdir> [--rvdb-taxmap <file>]

Expects under <outdir>: 05_genomad/<s>/<s>_virus_summary.tsv
                        06_bowtie2/<s>/<s>.coverage.tsv
                        07_diamond/<s>/{<s>.rvdb.tsv,<s>.uniref90.tsv}
                        09_viroid/<s>/<s>.viroid.tsv   (optional)
And under <db>/blastx/:  uniref90.taxmap.tsv
                         U-RVDBv32.0-prot.taxmap.tsv  (can override)
And optionally under <db>/blastn/: *_seqids.txt (viroid accession -> description)
EOF
}

sample= outdir= dbdir= rvdb_taxmap=
while [ $# -gt 0 ]; do
  case "$1" in
    -s) sample=$2; shift 2;;
    -o) outdir=$2; shift 2;;
    -d) dbdir=$2; shift 2;;
    --rvdb-taxmap) rvdb_taxmap=$2; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown arg: $1" >&2; usage >&2; exit 1;;
  esac
done
[ -n "$sample" ] && [ -n "$outdir" ] && [ -n "$dbdir" ] || { usage >&2; exit 1; }

blx="$dbdir/blastx"
bln="$dbdir/blastn"
rvdb_srt="$blx/U-RVDBv32.0-prot.taxmap.srt"
ur90_srt="$blx/uniref90.taxmap.srt"
seqids_srt="$bln/viroid_seqids.srt"
rvdb_taxmap=${rvdb_taxmap:-"$blx/U-RVDBv32.0-prot.taxmap.tsv"}
ur90_taxmap="$blx/uniref90.taxmap.tsv"

r="$outdir/07_diamond/$sample/$sample.rvdb.tsv"
u="$outdir/07_diamond/$sample/$sample.uniref90.tsv"
g="$outdir/05_genomad/$sample/${sample}_virus_summary.tsv"
[ -f "$g" ] || g="$outdir/05_genomad/$sample/$sample.virus_summary.tsv"
[ -f "$g" ] || g="$outdir/05_genomad/${sample}_genomad/${sample}_summary/${sample}_virus_summary.tsv"
cov="$outdir/06_bowtie2/$sample/$sample.coverage.tsv"
v="$outdir/09_viroid/$sample/$sample.viroid.tsv"
out="$outdir/evidence/$sample.evidence.tsv"

[ -f "$r" ] || { echo "missing: $r" >&2; exit 1; }
[ -f "$u" ] || { echo "missing: $u" >&2; exit 1; }
[ -f "$g" ] || g=/dev/null
[ -f "$cov" ] || cov=/dev/null
[ -f "$v" ] || v=/dev/null

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# ---- one-time caches ------------------------------------------------------
if [ ! -s "$rvdb_srt" ]; then
  [ -f "$rvdb_taxmap" ] || { echo "missing: $rvdb_taxmap" >&2; exit 1; }
  echo "building $rvdb_srt ..." >&2
  awk -F'\t' '$1!="rvdb_id" {print}' "$rvdb_taxmap" \
    | sort -t$'\t' -k1,1 \
    | awk -F'\t' '!seen[$1]++' > "$rvdb_srt"
fi
[ -s "$ur90_srt" ] || {
  [ -f "$ur90_taxmap" ] || { echo "missing: $ur90_taxmap" >&2; exit 1; }
  echo "building $ur90_srt ..." >&2
  awk -F'\t' '$1!="uniref90_id" {print}' "$ur90_taxmap" \
    | sort -t$'\t' -k1,1 > "$ur90_srt"
}
# viroid seqids: "ACCESSION rest-of-line" (single space after accession).
# Optional — absent file just means NA descriptions downstream.
if [ ! -s "$seqids_srt" ] && [ -d "$bln" ]; then
  found_seqids=0
  for f in "$bln"/*_seqids.txt; do
    [ -f "$f" ] || continue
    found_seqids=1
  done
  if [ "$found_seqids" -eq 1 ]; then
    mkdir -p "$bln"
    echo "building $seqids_srt ..." >&2
    for f in "$bln"/*_seqids.txt; do
      [ -f "$f" ] || continue
      awk '{ acc=$1; sub(/^[^ \t]+[ \t]+/,""); gsub(/\t/," "); print acc "\t" $0 }' "$f"
    done | sort -t$'\t' -k1,1 | awk -F'\t' '!seen[$1]++' > "$seqids_srt"
  fi
fi

# seqids cache may legitimately be absent (no *_seqids.txt); make it exist so
# the V-stream lookup below always has a loadable (possibly empty) file.
if [ ! -f "$seqids_srt" ]; then
  mkdir -p "$bln" 2>/dev/null || true
  : > "$seqids_srt" 2>/dev/null || seqids_srt=/dev/null
fi

# ---- rvdb hits: sseqid-keyed file -> protein FASTA header taxmap ------------
# joined: key contig sseqid prot ntacc prod pident aln_len mismatch gapopen
#         qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp |
#         prot ntacc organism product
awk -F'\t' '{
  split($2,a,"|");
  print $2 "\t" $1 "\t" $2 "\t" a[3] "\t" a[5] "\t" a[6] "\t" \
        $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" \
        $13 "\t" $14 "\t" $15 "\t" $16
}' "$r" \
  | sort -t$'\t' -k1,1 \
  | join -t$'\t' -a1 -1 1 -2 1 - "$rvdb_srt" > "$tmp/rvdb.joined"

# ---- top UniRef90 hit per contig -> taxmap join ---------------------------
# joined: id  contig sseqid pident aln_len mismatch gapopen qstart qend
#         sstart send evalue bitscore qlen slen qcovhsp scovhsp |
#         taxid organism repid description
awk -F'\t' '!seen[$1]++ {print $2"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' "$u" \
  | sort -t$'\t' -k1,1 \
  | join -t$'\t' -a1 -1 1 -2 1 - "$ur90_srt" > "$tmp/ur90.joined"

# ---- viroid blastn hits: one row per HSP, sseqid -> seqids description -----
# blastn sseqid is pipe-wrapped (gb|KY110721.1|) when the DB was built with
# -parse_seqids; fall back to the raw id otherwise. The small seqids cache is
# loaded whole, so hit-stream order is preserved (join would reorder it).
# V line: V contig sseqid acc pident aln_len mismatch gapopen qstart qend
#         sstart send evalue bitscore qlen slen qcovhsp scovhsp description
viroid_pre="$tmp/viroid.joined"
: > "$viroid_pre"
if [ "$v" != /dev/null ] && [ -s "$v" ]; then
  awk -F'\t' -v OFS='\t' -v srt="$seqids_srt" '
    FILENAME==srt { d[$1]=substr($0, index($0,"\t")+1); next }
    {
      if (index($2,"|")) { split($2,a,"|"); acc=a[2] } else acc=$2;
      desc = (acc in d) ? d[acc] : "NA";
      print "V", $1, $2, acc, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, desc
    }
  ' "$seqids_srt" "$v" > "$viroid_pre"
fi

# ---- assemble -------------------------------------------------------------
# The streams are not sorted here. Only the intermediate files above are
# sorted because POSIX join requires sorted keys.
{
  awk -F'\t' '$1 != "seq_name" {print "G\t" $0}' "$g"
  awk -F'\t' '!/^#/ {print "C\t" $0}' "$cov"
  sed 's/^/U\t/' "$tmp/ur90.joined"
  sed 's/^/R\t/' "$tmp/rvdb.joined"
  cat "$viroid_pre"
} | awk -F'\t' -v OFS='\t' '
  function add_contig(k) {
    if (!(k in present)) { present[k]=1; order[++norder]=k }
  }
  function emit(k,line,has,vline,vhas,    n,nz,i,y) {
    for (i=1;i<=9;i++) ga[i]="NA"
    for (i=1;i<=19;i++) rv[i]="NA"
    for (i=1;i<=19;i++) ua[i]="NA"
    for (i=1;i<=6;i++) ca[i]="NA"
    for (i=1;i<=17;i++) vi[i]="NA"

    n=split(k,p,"_")
    len=(k in gl) ? gl[k] : "NA"
    covv="NA"
    if (n >= 3 && p[n-1]=="cov") { if (len=="NA") len=p[n-2]; covv=p[n] }
    if (k in g) split(g[k],ga,"\t")
    if (k in c) split(c[k],ca,"\t")
    if (k in u) split(u[k],ua,"\t")

    if (has) {
      nz=split(line,z,FS)
      rv[1]=z[4]; rv[2]=z[5]; rv[3]=z[6]; rv[4]=z[7]; rv[5]=z[8]
      rv[6]=z[9]; rv[7]=z[10]; rv[8]=z[11]; rv[9]=z[12]; rv[10]=z[13]
      rv[11]=z[14]; rv[12]=z[15]; rv[13]=z[16]; rv[14]=z[17]
      rv[15]=z[18]; rv[16]=z[19]; rv[17]=z[20]; rv[18]=z[21]
      if (nz >= 25) {
        rv[2]=z[22]; rv[3]=z[23]; rv[4]=z[25]
        rv[19]=z[24]
      }
    }
    if (vhas) {
      nz=split(vline,y,FS)
      vi[1]=y[3]; vi[2]=y[4]
      for (i=4;i<=17;i++) vi[i]=y[i+1]
      vi[3]=(nz>=19) ? y[19] : "NA";
      for (i=20;i<=nz;i++) vi[3]=vi[3] FS y[i]
    }
    print k,len,covv,
      ga[1],ga[2],ga[3],ga[4],ga[5],ga[6],ga[7],ga[8],ga[9],
      rv[1],rv[2],rv[3],rv[4],rv[5],rv[6],rv[7],rv[8],rv[9],rv[10],rv[11],rv[12],rv[13],rv[14],rv[15],rv[16],rv[17],rv[18],rv[19],
      ua[1],ua[2],ua[3],ua[4],ua[5],ua[6],ua[7],ua[8],ua[9],ua[10],ua[11],ua[12],ua[13],ua[14],ua[15],ua[16],ua[17],ua[18],ua[19],
      ca[1],ca[2],ca[3],ca[4],ca[5],ca[6],
      vi[1],vi[2],vi[3],vi[4],vi[5],vi[6],vi[7],vi[8],vi[9],vi[10],vi[11],vi[12],vi[13],vi[14],vi[15],vi[16],vi[17]
  }
  $1=="G" { add_contig($2); gl[$2]=$3; g[$2]=$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; next }
  $1=="C" { add_contig($2); c[$2]=$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9; next }
  $1=="U" {
    for (i=0;i<15;i++) ua[1+i]=$(4+i);
    for (i=19;i<=22;i++) { ua[15+i-18] = (NF>=i) ? $i : "NA" }
    u[$3]=ua[1];
    for (i=2;i<=19;i++) u[$3]=u[$3]"\t"ua[i];
    add_contig($3)
    next
  }
  $1=="R" {
    add_contig($3); rv_count[$3]++; rv_line[$3,rv_count[$3]]=$0
  }
  $1=="V" {
    add_contig($2); vi_count[$2]++; vi_line[$2,vi_count[$2]]=$0
  }
  END {
    for (i=1;i<=norder;i++) {
      k=order[i]
      if (rv_count[k]) for (j=1;j<=rv_count[k];j++) emit(k,rv_line[k,j],1,"",0)
      if (vi_count[k]) for (j=1;j<=vi_count[k];j++) emit(k,"",0,vi_line[k,j],1)
      if (!rv_count[k] && !vi_count[k]) emit(k,"",0,"",0)
    }
  }
' > "$tmp/body"

mkdir -p "$(dirname "$out")"
cat > "$out" <<'HDR'
contig	contig_length	contig_cov	genomad_topology	genomad_coordinates	genomad_n_genes	genomad_genetic_code	genomad_virus_score	genomad_fdr	genomad_n_hallmarks	genomad_marker_enrichment	genomad_taxonomy	rvdb_sseqid	rvdb_protein_acc	rvdb_nt_acc	rvdb_product	rvdb_pident	rvdb_aln_len	rvdb_mismatch	rvdb_gapopen	rvdb_qstart	rvdb_qend	rvdb_sstart	rvdb_send	rvdb_evalue	rvdb_bitscore	rvdb_qlen	rvdb_slen	rvdb_qcovhsp	rvdb_scovhsp	rvdb_organism	uniref90_sseqid	uniref90_pident	uniref90_aln_len	uniref90_mismatch	uniref90_gapopen	uniref90_qstart	uniref90_qend	uniref90_sstart	uniref90_send	uniref90_evalue	uniref90_bitscore	uniref90_qlen	uniref90_slen	uniref90_qcovhsp	uniref90_scovhsp	uniref90_taxid	uniref90_organism	uniref90_repid	uniref90_description	bowtie2_numreads	bowtie2_covbases	bowtie2_coverage	bowtie2_meandepth	bowtie2_meanbaseq	bowtie2_meanmapq	viroid_sseqid	viroid_accession	viroid_description	viroid_pident	viroid_aln_len	viroid_mismatch	viroid_gapopen	viroid_qstart	viroid_qend	viroid_sstart	viroid_send	viroid_evalue	viroid_bitscore	viroid_qlen	viroid_slen	viroid_qcovhsp	viroid_scovhsp
HDR
cat "$tmp/body" >> "$out"
echo "wrote $out ($(($(wc -l < "$out")-1)) rows)"
