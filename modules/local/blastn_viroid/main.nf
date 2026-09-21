// blastn contigs against a pre-built viroid nucleotide database.
// Targets contigs in the 200–500 bp range where viroids (246–401 nt) will land.
process BLASTN_VIROID {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    path db_files

    output:
    tuple val(meta), path("${meta.id}.viroid.tsv"), emit: viroid

    script:
    """
    if [ -s ${contigs} ]; then
      blastn -db viroid_all_09_25_26_db -query ${contigs} -out ${meta.id}.viroid.tsv \\
        -task blastn -evalue 1e-5 -max_target_seqs 5 \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp' \\
        -num_threads ${task.cpus}
    else
      : > ${meta.id}.viroid.tsv
    fi
    """
}
