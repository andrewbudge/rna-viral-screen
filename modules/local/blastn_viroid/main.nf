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
      # 15 fields. blastn has no subject-coverage tag: unlike diamond it does not
      # accept scovhsp, and it silently DROPS an unrecognised trailing tag instead
      # of erroring, so requesting 16 fields yields 15 and the missing column only
      # shows up as empty evidence later. Do not re-add scovhsp here.
      blastn -db viroid_all_09_25_26_db -query ${contigs} -out ${meta.id}.viroid.tsv \\
        -task blastn -evalue 1e-5 -max_target_seqs 5 \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp' \\
        -num_threads ${task.cpus}
      # Guard: a short file means an outfmt tag was silently dropped. Must exit 0
      # when the file is empty (no hits is the common case), so no `&&` short-circuit.
      awk -F'\t' 'NF && NF!=15 {print FILENAME": line "NR" has "NF" fields, expected 15" > "/dev/stderr"; exit 1}' ${meta.id}.viroid.tsv
    else
      : > ${meta.id}.viroid.tsv
    fi
    """
}
