// Filter SPAdes contigs by length and mean coverage, and prefix every header
// with the sample id. Emit separate viroid-sized and standard screening sets.
process FILTER {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.viroid_contigs.fasta"), emit: contigs_small
    tuple val(meta), path("${meta.id}.contigs.fasta"), emit: contigs_large
    // Full filtered set (>=200 nt): the screening legs that are size-agnostic
    // (UniRef90 confirmation, map-back) run on this, one task per sample.
    tuple val(meta), path("${meta.id}.filtered_contigs.fasta"), emit: contigs_all

    script:
    """
    seqkit fx2tab ${contigs} \\
      | awk -F'\\t' -v minlen=${params.contig_min_length} -v mincov=${params.contig_min_cov} -v id='${meta.id}' '
          { split(\$1, a, "_"); if (a[4] >= minlen && a[6] >= mincov) print ">" id "_" substr(\$1, 6) "\\n" \$2 }' \\
      > ${meta.id}.filtered_contigs.fasta

    seqkit seq --min-len 200 --max-len 450 ${meta.id}.filtered_contigs.fasta \\
      > ${meta.id}.viroid_contigs.fasta

    seqkit seq --min-len 1000 ${meta.id}.filtered_contigs.fasta \\
      > ${meta.id}.contigs.fasta
    """
}
