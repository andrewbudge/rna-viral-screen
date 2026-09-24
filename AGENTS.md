# AGENTS.md — laney_plant_virus_pipeline (project notes for the coding agent)

Plant virus DETECTION pipeline for RNA-seq, Nextflow (nf-core-lite layout).
Owner: Andrew Budge. README.md is .gitignored until release; treat this file
as the live design record.

## Pipeline (current)

```
fastp ─► sortmerna ─► spades --rnaviral ─► filter (len>=200, cov)
  └──────────┴──────────► multiqc_report.html
         filter/contigs ─► diamond U-RVDB-prot · diamond UniRef90 · geNomad
                            plus bowtie2 read mapping · blastn viroid (optional)
```

- `modules/local/` — one process per file: fastp, sortmerna, sortmerna_index
  (removed; index now prebuilt), spades, filter, multiqc. blastn deferred;
  blastn_viroid optional for viroid detection (contigs 200-500 bp).
- Conf layout: `conf/base.config` (resources/retry), `conf/modules.config`
  (container/publishDir/ext.args), `conf/slurm.config` (CHPC executor).
- `main.nf` validates `--db` layout up front: sortmerna ref + index +
  genomad/genomad_db must exist before any process is submitted.

## Interface

```
nextflow run andrewbudge/laney_plant_virus_pipeline -profile slurm|local \
  --data_dir <root> --input samplesheet.tsv --outdir <dir> --db <dir>
```

- TSV columns: `sample  fastq_1  fastq_2`. Relative FASTQ paths resolve against
  `--data_dir`; absolute paths pass through. Relative `--db` also resolves
  against `--data_dir`.
- No `-resume` in the launchers (deliberate, documented). Pin `-r` to a tag for
  reproducibility.
- Repo is PUBLIC on GitHub; `nextflow run owner/repo` needs no token. `-latest`
  fetches latest main; run with `-r main` normally.

## DB layout (contract)

```
<db>/
  sortmerna/smr_v4.3_default_db.fasta   # MUST keep this exact basename
  sortmerna/index/                      # prebuilt index; name-encodes the ref basename
  blastn/U-RVDBv32.0.*                  # deferred/optional; makeblastdb -dbtype nucl -parse_seqids
  blastn/viroid_all_09_25_26_db.*       # viroid blastn DB; makeblastdb -dbtype nucl -parse_seqids
  blastx/                               # U-RVDBv32.0-prot.fasta, uniref90.fasta.gz (raw)
  blastx/uniref90.dmnd                  # Diamond UniRef90 database
  blastx/uniref90.taxmap.tsv            # header lookup: id, taxid, organism, repid, description
  blastx/uniref90.taxmap.srt            # one-time cache (agg script sorts it); safe to delete
  blastx/U-RVDBv32.0-prot.taxmap.tsv    # protein FASTA header lookup: rvdb id, protein/nt accessions, organism, product
  blastx/U-RVDBv32.0-prot.taxmap.srt    # one-time sorted cache; safe to delete
  blastn/*_seqids.txt                   # viroid accession -> description ("ACC rest-of-line"); required by evidence merge
  genomad/genomad_db/                   # genomad download-database output
```

Index must match reference basename (Nextflow stages symlinked refs under the
link name → SortMeRNA index mismatch). Do NOT symlink the reference.

## Pinned images (AVX constraints — lonepeak is Sandy Bridge, AVX only)

| tool | image | why |
|---|---|---|
| fastp | fastp:0.23.4--h125f33a_5 | |
| sortmerna | sortmerna:4.3.4--h9ee0642_0 | 4.3.6+/7.0.0 need x86-64-v3 (AVX2+BMI2+FMA) → SIGILL |
| spades | spades:4.0.0--haf24da9_4 | 4.2+ needs AVX2 → SIGILL; 4.0.0 supports --rnaviral |
| seqkit | seqkit:2.13.0--he881be0_0 | Go binary, AVX-safe |
| blast | blast:2.17.0--hb02a186_1 | blastn deferred; 2.17.0 verified working on lonepeak |
| bowtie2 | bowtie2:2.5.5--ha27dd3b_0 | SSE-era code, AVX-safe |
| genomad | genomad:1.12.0--pyhdfd78af_0 | SIGILL risk: bundled mmseqs2/TF need modern CPUs; verify before production, may need non-lonepeak partition |
| multiqc | multiqc:1.25.2--pyhdfd78af_0 | |

Test any new container with `apptainer exec <img> <tool> --version` before
trusting it on lonepeak. tagfind (`./bin/tagfind <tool> [prefix]`) lists
biocontainer tags. diamond:2.2.6 flagged AVX2-risk — verify before use.

## Screening strategy (V1 LOCKED with user)

Goal: **curate evidence for experts to sift** — no prefiltration, accrue info.
Host-agnostic; can't guarantee completeness. Four active parallel legs on
FILTER's contigs, each feeds one raw per-leg TSV; blastn is deferred.

| leg | tool → DB | what it answers |
|---|---|---|
| 1 | blastn → U-RVDB (nucleotide) | (deferred) known-virus relative? (homology) |
| 2 | diamond blastx → U-RVDB-prot | divergent-virus homology (protein) (implemented) |
| 3 | diamond blastx → UniRef90 | unbiased: is the best hit viral? (confirmation, implemented) |
| 4 | geNomad | viral by sequence/marker signal, no homology needed (implemented, confirmed) |
| 5 | blastn → viroid DB | viroid detection (contigs 200-450 bp) (implemented) |
| support | bowtie2 map-back reads → contigs | per-contig read support and confidence for calls |

Decisions made:
- **UniRef90 only — no nr.** UniRef90 ⊂ nr for any known virus; nr is ~500M seqs
  (hours to format, huge disk) and adds nothing except for viruses absent from
  UniRef90. Confirmation leg = contigs themselves vs UniRef90, NOT the viral
  hit sequences.
- **geNomad only for V1 — no VirSorter2.** geNomad was retrained for RNA +
  giant viruses (99.9% of RdRP contigs flagged vs ~44% for others) and beats
  all tools in its benchmark (MCC 95.3%); handles <3kb contigs; ~26x faster
  than VirSorter2. Skip k-mer CNN tools (PPR-Meta/DeepVirFinder/Seeker) — they
  false-positive on eukaryotic seqs at 2-3x rate. VirSorter2 (`--include-groups
   RNA`) is the agreed backup if geNomad underperforms.
- **blastn is deferred while geNomad + bowtie2 map-back are the focus.**
- **Viroid detection via blastn** — RVDB/UniRef90/geNomad all miss viroids
  (non-coding RNA, 246-401 nt). Always-on leg using `<db>/blastn/viroid_all_09_25_26_db.*`
  (13,132 sequences from NCBI GenBank, 2026-09-15).
- Evidence aggregation implemented (`bin/aggregate_evidence.sh`): per-contig ×
  RVDB-hit join of geNomad + diamond RVDB/UniRef90 + protein header taxmaps + bowtie2.
- **Map-back support is a confidence leg.** rRNA-depleted reads mapped back to
  filtered contigs support calls such as the RNA mitovirus; BAM enables
  `samtools depth`/`idxstats` for the deferred evidence table.
- **Screening-leg channel rule (v0.2.3 fix):** FILTER emits three sets —
  `contigs_all` (>=200nt, the union) drives UniRef90 + bowtie2 map-back;
  `contigs_large` (>=1000nt) drives RVDB blastx + geNomad; `contigs_small`
  (200-450bp) drives viroid blastn. Previously UniRef90/bowtie2 were fed a
  `mix(small, large)` of two emissions per sample → two tasks writing the same
  output name, publishDir overwrite, and duplicate join keys in EVIDENCE.
- **Evidence table fixes in v0.2.3:** bowtie2 columns were off by one (took
  `covbases..` where samtools coverage puts `numreads`), and the table grew
  56→73 cols with the viroid leg (one row per HSP, sseqid + accession +
  seqids.txt description + 14 blast fields).
- **GENOMAD takes `val db`, not `path db` (2026-09-24).** `path db` stages a
  symlink into the work dir and renders `${db}` *relative*, so any `rm` or
  `find -type f` aimed at `genomad_db/` from inside a task walks straight
  through into the shared database. The 2026-09-22 incident left
  `<db>/genomad/genomad_db` holding exactly the 8 `genomad_mini_db*` symlinks
  the tarball ships and zero regular files — a `-type f` delete signature,
  not a staging bug. `val db` stages nothing and renders the absolute path.
  Two guards now catch the two failure modes separately: the task asserts
  `${db}/version.txt` *inside the container* (unbound path → message in
  `<sample>.genomad.log`), and `main.nf` asserts `version.txt` + `genomad_db`
  + `nodes.dmp` + `names.dmp` + `genomad_marker_metadata.tsv` at launch
  (deleted DB → fails before any SLURM job). `conf/slurm.config`'s
  `runOptions` bind now reaches geNomad on its own — do not drop it.
  DIAMOND_* stay on `path db` (single `.dmnd` files, harmless) and SORTMERNA
  stays on `path ref` (the index-basename constraint above depends on the
  staged link name).

DBs needed (sources already staged in `<db>/blastx/` + user downloading):
- U-RVDB v32.0 (rvdb.dbi.udel.edu), deferred → `makeblastdb -dbtype nucl -parse_seqids`
- U-RVDB-prot (Institut Pasteur, rvdb-prot.pasteur.fr) → `diamond makedb --in U-RVDBv32.0-prot.fasta -d U-RVDB-prot`
- UniRef90 (uniref90.fasta.gz on HPC) → `diamond makedb`
- geNomad DB → `genomad download-database` on HPC

Note: RVDB includes endogenous retrovirus / LTR-retrotransposon sequences —
plant hosts will hit them. That's expected; UniRef90 best-hit sorts it.
RVDB "Putative Non Viral Annotation" file exists if we want to tag hits later.

Notes: blastx >> blastn for sensitivity on divergent viruses (protein diverges
~10x slower). "Both legs viral" = confidence flag, NOT a hard filter. Diamond =
100-1000x faster than blastx. Diamond 2.2.6 AVX2-risk on lonepeak — verify
`diamond --version` before use.

## Output

```
<outdir>/
  01_fastp/<sample>/  02_sortmerna/<sample>/  03_spades/<sample>/
  04_filter/<sample>/{<sample>.filtered_contigs.fasta,<sample>.viroid_contigs.fasta,<sample>.contigs.fasta}
                                               # filtered: >=200nt cov>=10 (screening set); viroid: 200-450 bp; standard: >=1000 bp
  05_genomad/<sample>/{<sample>_virus_summary.tsv,<sample>_virus.fna,<sample>_virus_proteins.faa,<sample>_summary.json,*.genomad.log}
  06_bowtie2/<sample>/{*.sorted.bam,*.sorted.bam.bai,*.coverage.tsv,*.bowtie2.log}
  07_diamond/<sample>/{<sample>.rvdb.tsv,<sample>.uniref90.tsv}
  08_multiqc/  evidence/<sample>.evidence.tsv   # 73-col joined table via bin/aggregate_evidence.sh
  09_viroid/<sample>/<sample>.viroid.tsv    # blastn vs viroid DB (contigs 200-450 bp), merged into evidence
  pipeline_info/{timeline,report,trace,dag}
```

Evidence aggregation: `bin/aggregate_evidence.sh -s <sample> -o <outdir> -d <dbdir>`
joins geNomad + diamond RVDB/UniRef90 + protein header taxmaps + bowtie2 coverage +
viroid blastn on contig name into a 73-col TSV (one row per union contig; repeated per
RVDB hit and per viroid HSP when multiple hits exist; each multi-row carries NA for the
other multi-row leg). One-time sorted caches (`.srt`) auto-built in `<db>/blastx/` and
`<db>/blastn/` on first run. Final rows retain stream order; only join inputs and lookup
caches are sorted as required by `join`. The viroid leg joins sseqid (pipe-wrapped
`gb|ACC|`) to `<db>/blastn/*_seqids.txt` for `viroid_accession` + `viroid_description`.

Diamond outfmt carries only the subject seqid — the first header token before a
space. RVDB-prot descriptions thus truncate to one word (e.g. `essential`,
`Gag-Pol`). So each diamond leg gets a **header lookup TSV** (one-time build,
stored beside the DB) for the evidence aggregation to join on sseqid:

| leg | lookup key | file |
|---|---|---|
| UniRef90 | `UniRef90_*` id | `<db>/blastx/uniref90.taxmap.tsv` (id, taxid, organism, repid, description) — built |
| U-RVDB-prot | complete DIAMOND `sseqid` | `<db>/blastx/U-RVDBv32.0-prot.taxmap.tsv` (rvdb id, protein accession, nucleotide accession, organism, product) — built from protein FASTA headers |

Both DIAMOND legs use a custom 16-field outfmt: the standard 12 fields plus
`qlen`, `slen`, `qcovhsp`, and `scovhsp`.

RVDB additionally ships `RVDB_AnnotationList_Current.tab.gz` (non-viral/ERV
tags) if we want to tag Gag-Pol/LTR hits in aggregation later. RVDB protein
headers do not carry taxids or lineages, so those columns are omitted until a
reliable protein-accession taxonomy mapping is available.

blastn remains deferred; its planned outfmt is deliberately 12 columns.
qcovs = fraction of contig matched, `length` vs `slen` = does the contig span
the whole genome. staxids/sscinames stay blank without `makeblastdb -taxid_map`.

## Potential add-ons after 0.1.0

- geNomad `annotate/*_genes.tsv` aggregation: marker count, USCG count,
  plasmid/virus hallmark counts, and hallmark identities, coordinates, scores,
  taxonomy, and descriptions for every contig. Use the all-contig annotation
  file rather than only `summary/*_virus_genes.tsv` to preserve no-prefilter
  evidence collection. Do not interpret marker `taxname` alone as viral.
- Depth uniformity from `samtools depth` or `mosdepth`: median depth, depth CV,
  uncovered bases, and fractions at >=10x and >=100x.
- Genome completeness evidence: ORFfinder/EMBOSS getorf, HMMER/Pfam or
  InterProScan domains, reference gene complement, Bandage assembly graphs,
  and experimental RACE for RNA termini. CheckV is mainly applicable to DNA
  viruses and phages, not most RNA viruses.
- Novelty and taxonomy: PalmScan/PalmDB or RdRP profile HMMs, MMseqs2 protein
  clustering, MAFFT alignments, IQ-TREE 2 phylogenies, and ICTV family-specific
  demarcation criteria.
- Host association and active infection: sample metadata, multi-sample
  co-abundance, broad Kraken2/Kaiju/DIAMOND taxonomy, strand-specific RNA-seq,
  21-24 nt viral small-RNA signatures, and RT-qPCR/ddPCR or in situ validation.
- Contamination controls: extraction blanks, technical/biological replicates,
  cross-sample prevalence and index-hopping checks, and `decontam` analysis.
- Absence confidence: positive spike-ins, defined detection limits, read
  downsampling, translated/profile-HMM sensitivity benchmarks, and replicate
  consistency. A negative screen cannot establish true virus absence.
- **[PARKED 2026-09-22, pending PI feedback]** General-purpose viral genome
  annotation CLI (Rust, clap, static): reference-transfer annotation +
  readthrough-aware ORF calling (126K/183K class) + plant genetic code + GFF3/
  .tbl/table2asn emission + MIUViG reporting fields. Gap evidence in
  `~/literature/viral-annotation/` (VAPiD: Shean 2019 PMID 30674273 —
  Python-2, human-virus-only, CDS-only, readthrough untested; Prodigal: Hyatt
  2010 PMID 20211023 — start→first-stop model cannot emit readthrough CDS;
  practice anchors: Zisi 2024 PMID 38774217, Esmaeilzadeh 2023 PMID 37894095).

## Verification workflow

Full local runtime tests with synthetic data under /tmp (proven pattern):
random 3kb genome → sortmerna ref + index; second random 3kb genome → 500 FR
151bp reads (p in [0,L-302], r1=B[p:p+151], r2=rc(B[p+151:p+302])) + assembly
target. Validate publishDir closures + output declarations only at runtime
(preview doesn't catch them). Always: `nextflow run -preview`, `git diff --check`.

## HPC

- lonepeak SLURM, account ogden. Images cached in $PVP/apptainer_nf via
  NXF_APPTAINER_CACHEDIR. PVP=/scratch/general/vast/abudge/plant_virus_pipeline_storage.
- Work dir $PVP/nf_work, outdir $PVP/qc_assembly. Production run:
  `nextflow run andrewbudge/laney_plant_virus_pipeline -latest -r main -resume
  -profile slurm -work-dir "$PVP/nf_work" --data_dir "$PVP"
  --input "$HOME/samplesheet.tsv" --outdir "$PVP/qc_assembly" --db databases`
- Resources measured on Laney18 (18.2M pairs): FASTP 2m/3GB, SORTMERNA 1h24m/4.5GB.
  SPADES seeds: 8 cpus/32GB/12h. Log real runtimes in base.config comments.
