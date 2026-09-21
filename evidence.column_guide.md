# Evidence Table — Column Guide

Companion to the pipeline evidence table (`<outdir>/evidence/<sample>.evidence.tsv`,
73 columns). One row per contig, repeated per U-RVDB-prot hit and per viroid
blastn HSP when those hits exist (geNomad, UniRef90 and Bowtie2 values repeat
across rows of the same contig). A contig with both RVDB and viroid hits gets
RVDB rows and viroid rows; each row is `NA` for the other leg (in practice the
legs target disjoint size classes). `NA` means no data / no hit for that leg.
Example values below are illustrative, taken from a mitovirus candidate contig.

How to read it: **RVDB** = divergent-virus homology leg (DIAMOND blastx vs
U-RVDB-prot). **UniRef90** = unbiased confirmation leg (DIAMOND blastx, best hit
only). **geNomad** = homology-free viral call. **Bowtie2** = read support for
the contig. **Viroid** = blastn vs the viroid nucleotide DB (200–450 bp
contigs). A contig where both DIAMOND legs agree on viral is the strongest
candidate; a high-identity full-length Viroid row is a viroid call.

## Contig (3)

| Column | What it means | Example |
|---|---|---|
| `contig` | Assembled contig name; encodes sample, ID, length and k-mer coverage. | `<sample>_<id>_length_<len>_cov_<cov>` |
| `contig_length` | Contig length in nucleotides. | 2770 |
| `contig_cov` | Assembler k-mer coverage (assembly depth proxy). | 1934.4 |

## geNomad — homology-free viral call (9)

| Column | What it means | Example |
|---|---|---|
| `genomad_topology` | Provirus topology call (e.g. No terminal repeats). | No terminal repeats |
| `genomad_coordinates` | Provirus coordinates if a provirus region called, else NA. | NA |
| `genomad_n_genes` | Number of predicted genes on the contig. | 2 |
| `genomad_genetic_code` | Genetic code used for gene prediction. | 11 |
| `genomad_virus_score` | geNomad virus score, 0–1; higher = more viral. | 0.9993 |
| `genomad_fdr` | False-discovery rate estimate; lower = more confident. | 0.0006 |
| `genomad_n_hallmarks` | Count of viral hallmark markers found. | 1 |
| `genomad_marker_enrichment` | Viral marker enrichment score. | 1.7183 |
| `genomad_taxonomy` | geNomad taxonomic assignment, semicolon-delimited ranks. | …Mitoviridae… |

## RVDB — DIAMOND blastx vs U-RVDB-prot, up to 5 hits per contig (19)

| Column | What it means | Example |
|---|---|---|
| `rvdb_sseqid` | Raw DIAMOND subject ID (full pipe string); provenance / join key. | `acc\|GENBANK\|UJQ92534.1\|…\|putative` |
| `rvdb_protein_acc` | Protein accession parsed from sseqid / taxmap. | UJQ92534.1 |
| `rvdb_nt_acc` | Source nucleotide accession. | MZ679949 |
| `rvdb_product` | Protein product; taxmap text is fuller than the sseqid tail. | putative RNA-dependent RNA polymerase |
| `rvdb_pident` | Percent identity of the alignment. | 40.7 |
| `rvdb_aln_len` | Alignment length in amino acids. | 734 |
| `rvdb_mismatch` | Mismatches in the alignment. | 410 |
| `rvdb_gapopen` | Gap openings in the alignment. | 13 |
| `rvdb_qstart` / `rvdb_qend` | Hit coordinates on contig (translated frame); start > end = minus strand. | 2462 / 273 |
| `rvdb_sstart` / `rvdb_send` | Hit coordinates on the subject protein. | 14 / 726 |
| `rvdb_evalue` | Expect value; smaller = more significant. | 6.92e-171 |
| `rvdb_bitscore` | Bit score; larger = better hit. | 524 |
| `rvdb_qlen` | Query (contig) length in nt. | 2770 |
| `rvdb_slen` | Subject protein length in aa. | 754 |
| `rvdb_qcovhsp` | % of contig covered by the hit (HSP). | 79.1 |
| `rvdb_scovhsp` | % of subject protein covered by the hit. | 94.6 |
| `rvdb_organism` | Subject organism from RVDB taxmap (no taxid/lineage available). | Mitoviridae sp. |

## UniRef90 — DIAMOND blastx, best hit only (19)

| Column | What it means | Example |
|---|---|---|
| `uniref90_sseqid` | UniRef90 cluster ID of best hit. | UniRef90_A0ABM9WIR1 |
| `uniref90_pident` | Percent identity. | 37.3 |
| `uniref90_aln_len` | Alignment length in aa. | 742 |
| `uniref90_mismatch` | Mismatches. | 436 |
| `uniref90_gapopen` | Gap openings. | 12 |
| `uniref90_qstart` / `uniref90_qend` | Hit coordinates on contig; start > end = minus strand. | 2435 / 273 |
| `uniref90_sstart` / `uniref90_send` | Hit coordinates on subject. | 7 / 740 |
| `uniref90_evalue` | Expect value. | 6.35e-144 |
| `uniref90_bitscore` | Bit score. | 455 |
| `uniref90_qlen` | Query (contig) length in nt. | 2770 |
| `uniref90_slen` | Subject length in aa. | 782 |
| `uniref90_qcovhsp` | % of contig covered. | 78.1 |
| `uniref90_scovhsp` | % of subject covered. | 93.9 |
| `uniref90_taxid` | NCBI taxid of the UniRef90 entry. | 2080459 |
| `uniref90_organism` | Organism of best hit — host vs virus sorts true/false positives. | Dahlia pinnata mitovirus 1 |
| `uniref90_repid` | Representative member ID of the cluster. | A0ABM9WIR1_9VIRU |
| `uniref90_description` | Protein description. | RNA-dependent RNA polymerase |

## Bowtie2 — read map-back support (6)

| Column | What it means | Example |
|---|---|---|
| `bowtie2_numreads` | Reads mapped back to the contig. | 95304 |
| `bowtie2_covbases` | Contig bases covered by ≥1 read. | 2770 |
| `bowtie2_coverage` | Breadth of coverage in %; 100 = fully covered. | 100 |
| `bowtie2_meandepth` | Mean read depth across the contig. | 4645.14 |
| `bowtie2_meanbaseq` | Mean base quality of mapped reads. | 38.9 |
| `bowtie2_meanmapq` | Mean mapping quality. | 41.8 |

## Viroid — blastn vs the viroid nucleotide DB, one row per HSP (17)

Filled only on rows from the viroid leg; RVDB rows and plain contig rows carry
NA here, and viroid rows carry NA in the RVDB columns.

| Column | What it means | Example |
|---|---|---|
| `viroid_sseqid` | Raw blastn subject ID (pipe-wrapped with `-parse_seqids`); join key / provenance. | `gb\|KY110721.1\|` |
| `viroid_accession` | Accession parsed from the sseqid. | KY110721.1 |
| `viroid_description` | Description from `<db>/blastn/*_seqids.txt`. | ...viroid complete genome |
| `viroid_pident` | Nucleotide percent identity. | 96.0 |
| `viroid_aln_len` | Alignment length in nt. | 375 |
| `viroid_mismatch` | Mismatches in the alignment. | 10 |
| `viroid_gapopen` | Gap openings. | 4 |
| `viroid_qstart` / `viroid_qend` | Hit coordinates on contig; start > end = minus strand. | 55 / 427 |
| `viroid_sstart` / `viroid_send` | Hit coordinates on the subject viroid genome. | 372 / 1 |
| `viroid_evalue` | Expect value. | 5.04e-171 |
| `viroid_bitscore` | Bit score; larger = better hit. | 596 |
| `viroid_qlen` | Query (contig) length in nt. | 446 |
| `viroid_slen` | Subject viroid genome length in nt (246–401 typical). | 372 |
| `viroid_qcovhsp` | % of contig covered by the HSP. | 84.0 |
| `viroid_scovhsp` | % of subject viroid covered by the HSP. | 99.0 |

A near-full-length, high-identity pair (`qcovhsp` and `scovhsp` both high,
`viroid_slen` matching `contig_length`) is the classic complete-viroid call:
viroids are 246–401 nt circular non-coding RNAs, so the contig should span the
whole subject. Partial or multiple HSPs per contig can mean a variant, a
truncated assembly, or a circular permutation; check `qstart/qend` vs
`sstart/send`.

## Notes

- `rvdb_protein_acc` / `rvdb_nt_acc` are pipe-fields of `rvdb_sseqid`;
  `rvdb_product` and `rvdb_organism` come only from the header taxmap.
  Keep `rvdb_sseqid` in the full table as the join key / provenance that traces
  a row back to the raw DIAMOND output and the source FASTA header; hide it in
  review views.
- UniRef90 keeps the single best hit per contig; RVDB keeps up to 5 hits meeting
  the query-cover threshold. One contig can therefore occupy several rows that
  differ only in the RVDB columns.
- `qcovhsp` = fraction of the contig matched; `aln_len` vs `slen` shows whether
  the contig spans the whole subject. A short high-identity hit with low `qcovhsp`
  is a domain-level match, not a whole-gene call.
- Strand: `qstart > qend` means the hit is on the minus strand (translated frame
  on the reverse complement). This is normal and carries no quality meaning.

### Worked example 1 — strong candidate (both legs viral)

Contig `Laney18_1313` (2770 nt, cov 1934.4):

- geNomad: score 0.9993, FDR 0.0006, taxonomy …Mitoviridae… — viral call with a
  family assignment.
- RVDB: 40.7% identity over 734 aa to `UJQ92534.1` (putative RdRp, Mitoviridae
  sp.), e-value 6.92e-171, qcov 79.1 / scov 94.6 — a near-full-protein remote
  homolog, exactly the signal the RVDB leg exists to catch.
- UniRef90: best hit is itself viral — 37.3% to Dahlia pinnata mitovirus 1 RdRp
  (e-value 6.35e-144), not a plant protein. Both legs agree.
- Bowtie2: 95,304 reads, 100% breadth, depth 4645x — deep, even support, not a
  low-coverage artifact.

Verdict: prioritize. Divergent mitovirus; worth phylogenetic follow-up.

### Worked example 2 — host false positive (legs disagree)

Contig `Laney18_4083` (1897 nt, cov 19.1):

- geNomad: score 0.9990 yet taxonomy Unclassified, 1 gene, 1 hallmark, enrichment
  1.72 — a high score with nothing behind it.
- RVDB: no hit at all.
- UniRef90: 98.6% identity over the full 559 aa protein to `B9HUU6_POPTR`, a
  poplar MFS transporter (e-value 0.0) — a near-identical host gene.
- Bowtie2: 642 reads, 100% breadth, depth 45x — genuinely expressed, which is
  why it assembled, but expression is not viral evidence.

Verdict: deprioritize. Single hallmark firing on a host transporter; the UniRef90
leg overrules geNomad. This disagreement is the design working — the conflict
stays visible in the table instead of being silently filtered.

### Worked example 3 — viroid call (from output_v0.2.2 batch 1)

Contig `Laney29_47263_length_446_cov_53966.59` (446 nt, k-mer cov ~54k):

- Viroid leg: best hit `gb|KY110721.1|` — 96.0% identity over 375 nt
  (e-value 5.04e-171), a second HSP covers contig 1–54; other subjects
  (`KY110722.1`, `PV339791.1`, ...) hit at 95–97% over ~340–372 nt.
- Subject genomes are ~370 nt (`viroid_slen`) against a 446 nt contig
  (`viroid_qlen`), near-full `scovhsp` — classic viroid size class; the contig
  is longer than the genome because circular viroids assemble into a
  circularly-permuted concatemer, hence the split/tiled HSPs.
- geNomad / RVDB columns NA: the contig never entered those legs (200–450 bp
  and non-coding), and UniRef90 should agree with NA — viroids encode no
  protein.
- Depth: k-mer cov 54k = extremely high; viroids replicate hard.

Verdict: prioritize as a viroid candidate. Follow up by circularizing the
contig in silico and re-blasting, then compare against the reference variant
collection for isolate-level calls.
