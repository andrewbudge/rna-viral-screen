// Merge all evidence legs into one descriptive per-sample TSV.
process EVIDENCE {
    tag "${meta.id}"
    label 'assembly'

    input:
    tuple val(meta), path(rvdb), path(uniref90), path(virus_summary), path(coverage), path(viroid)
    path rvdb_taxmap
    path uniref90_taxmap
    path viroid_seqids
    path aggregate_script

    output:
    tuple val(meta), path("${meta.id}.evidence.tsv"), emit: evidence

    script:
    """
    mkdir -p inputs/05_genomad/${meta.id} inputs/06_bowtie2/${meta.id} inputs/07_diamond/${meta.id} inputs/09_viroid/${meta.id} db/blastx db/blastn
    cp ${virus_summary} inputs/05_genomad/${meta.id}/${meta.id}_virus_summary.tsv
    cp ${coverage} inputs/06_bowtie2/${meta.id}/${meta.id}.coverage.tsv
    cp ${rvdb} inputs/07_diamond/${meta.id}/${meta.id}.rvdb.tsv
    cp ${uniref90} inputs/07_diamond/${meta.id}/${meta.id}.uniref90.tsv
    cp ${viroid} inputs/09_viroid/${meta.id}/${meta.id}.viroid.tsv
    ln -s ../../${rvdb_taxmap} db/blastx/U-RVDBv32.0-prot.taxmap.tsv
    ln -s ../../${uniref90_taxmap} db/blastx/uniref90.taxmap.tsv
    for f in ${viroid_seqids}; do ln -s ../../\$f db/blastn/\$(basename \$f); done
    bash ${aggregate_script} -s ${meta.id} -o inputs -d db
    mv inputs/evidence/${meta.id}.evidence.tsv ${meta.id}.evidence.tsv
    """
}
