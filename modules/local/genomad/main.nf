// Classify filtered contigs with geNomad. Plasmid_* outputs also exist but are
// not declared: V1 is virus-focused only.
process GENOMAD {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(contigs)
    val db

    output:
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus_summary.tsv"), emit: virus_summary
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus.fna"), emit: virus_contigs
    tuple val(meta), path("${meta.id}_genomad/${meta.id}_summary/${meta.id}_virus_proteins.faa"), emit: virus_proteins
    path "${meta.id}_genomad/${meta.id}_summary/${meta.id}_summary.json", emit: summary_json
    path "${meta.id}.genomad.log", emit: log

    script:
    """
    if [ ! -f ${db}/version.txt ]; then
        echo "geNomad database not reachable inside container: ${db}/version.txt" | tee /dev/stderr > ${meta.id}.genomad.log
        exit 1
    fi
    cp ${contigs} ${meta.id}.fna
    genomad end-to-end --cleanup --threads ${task.cpus} --sensitivity 5.0 --lenient-taxonomy --full-ictv-lineage --enable-score-calibration ${meta.id}.fna ${meta.id}_genomad ${db} 2> ${meta.id}.genomad.log
    """
}
