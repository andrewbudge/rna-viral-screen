#!/usr/bin/env nextflow
//
// plant-virus QC + assembly + viral screening: adapter/quality trim (fastp),
// rRNA depletion (SortMeRNA), de novo RNA viral assembly (SPAdes --rnaviral),
// contig filtering, DIAMOND U-RVDB-prot and UniRef90, geNomad, bowtie2 read
// mapping, plus a MultiQC report per run.
//
// Processes live in modules/local/; their containers, publishDir and tool
// arguments are configured in conf/modules.config. Resources/retry live in
// conf/base.config, the HPC executor in conf/slurm.config.
//
//   nextflow run main.nf -profile local --input samplesheet.tsv --data_dir <dir> --db <dir>
//   nextflow run main.nf -profile slurm --input samplesheet.tsv --data_dir <dir> --db <dir>
//
include { FASTP          } from './modules/local/fastp/main.nf'
include { SORTMERNA      } from './modules/local/sortmerna/main.nf'
include { SPADES         } from './modules/local/spades/main.nf'
include { FILTER         } from './modules/local/filter/main.nf'
include { BOWTIE2        } from './modules/local/bowtie2/main.nf'
include { DIAMOND_RVDB   } from './modules/local/diamond_rvdb/main.nf'
include { DIAMOND_UNIREF90 } from './modules/local/diamond_uniref90/main.nf'
include { GENOMAD        } from './modules/local/genomad/main.nf'
include { MULTIQC        } from './modules/local/multiqc/main.nf'
include { SAMTOOLS_SORT  } from './modules/local/samtools_sort/main.nf'
include { EVIDENCE       } from './modules/local/evidence/main.nf'
include { BLASTN_VIROID  } from './modules/local/blastn_viroid/main.nf'

// Relative paths resolve against a caller-supplied root; absolute paths pass
// through untouched.
def resolve(path, root, what) {
    if (!path) error "${what} is not set"
    def p = path.toString().trim()
    if (p.startsWith('/')) return file(p, checkIfExists: true)
    if (!root) error "${what} is relative ('${p}') but its root is not set"
    return file("${root.toString().replaceAll('/$', '')}/${p}", checkIfExists: true)
}

workflow {
    if (params.data_dir && !file(params.data_dir).exists())
        error "--data_dir does not exist: ${params.data_dir}"

    db = resolve(params.db, params.data_dir, '--db')
    if (!db.isDirectory()) error "--db is not a directory: ${db}"

    sortmerna_ref = file("${db}/sortmerna/smr_v4.3_default_db.fasta", checkIfExists: true)
    sortmerna_idx = file("${db}/sortmerna/index", checkIfExists: true)
    if (!sortmerna_idx.isDirectory())
        error "SortMeRNA index is not a directory: ${sortmerna_idx}"

    genomad_db = file("${db}/genomad/genomad_db")
    if (!genomad_db.isDirectory())
        error "geNomad database not found: ${genomad_db} — run 'genomad download-database <db>/genomad' to create it"
    genomad_required = ['version.txt', 'genomad_db', 'nodes.dmp', 'names.dmp', 'genomad_marker_metadata.tsv']
    genomad_missing = genomad_required.findAll { !new File(genomad_db.toString(), it).exists() }
    if (genomad_missing)
        error "geNomad database incomplete at ${genomad_db}: missing ${genomad_missing.join(', ')} — run 'genomad download-database <db>/genomad' to rebuild it"

    rvdb_prot_db = file("${db}/blastx/U-RVDB-prot.dmnd")
    uniref90_db = file("${db}/blastx/uniref90.dmnd")
    rvdb_taxmap = file("${db}/blastx/U-RVDBv32.0-prot.taxmap.tsv", checkIfExists: true)
    uniref90_taxmap = file("${db}/blastx/uniref90.taxmap.tsv", checkIfExists: true)
    viroid_db_files = file("${db}/blastn/viroid_all_09_25_26_db.*", checkIfExists: true)
    if (!viroid_db_files.find { it.name.endsWith('.nin') } || !viroid_db_files.find { it.name.endsWith('.nsq') })
        error "viroid BLAST database is incomplete: ${db}/blastn/viroid_all_09_25_26_db.*"
    viroid_seqids = file("${db}/blastn/*_seqids.txt")
    if (!viroid_seqids)
        error "viroid accession->description lookup missing: ${db}/blastn/*_seqids.txt"
    aggregate_script = file("${baseDir}/bin/aggregate_evidence.sh", checkIfExists: true)

    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "samplesheet ${params.input}: every row needs sample, fastq_1, fastq_2 — got: ${row}"
            }
            tuple([id: row.sample],
                  [resolve(row.fastq_1, params.data_dir, "${row.sample}.fastq_1"),
                   resolve(row.fastq_2, params.data_dir, "${row.sample}.fastq_2")])
        }
        .toList()
        .map { rows ->
            if (!rows) error "samplesheet ${params.input} has no data rows"
            // Duplicates would collide in publishDir and overwrite each other.
            def dupes = rows.collect { it[0].id }.countBy { it }.findAll { it.value > 1 }.keySet()
            if (dupes) error "duplicate sample names in ${params.input}: ${dupes}"
            rows
        }
        .flatMap()

    ch_ref = Channel.value(sortmerna_ref)
    ch_idx = Channel.value(sortmerna_idx)

    FASTP(ch_samples)
    SORTMERNA(FASTP.out.reads, ch_ref, ch_idx)
    SPADES(SORTMERNA.out.clean)
    FILTER(SPADES.out.contigs)
    DIAMOND_RVDB(FILTER.out.contigs_large, Channel.value(rvdb_prot_db))
    DIAMOND_UNIREF90(FILTER.out.contigs_all, Channel.value(uniref90_db))
    GENOMAD(FILTER.out.contigs_large, Channel.value(genomad_db))
    BOWTIE2(SORTMERNA.out.clean.join(FILTER.out.contigs_all))
    BLASTN_VIROID(FILTER.out.contigs_small, Channel.value(viroid_db_files))
    SAMTOOLS_SORT(BOWTIE2.out.sam)
    evidence_inputs = DIAMOND_RVDB.out.rvdb
        .join(DIAMOND_UNIREF90.out.uniref90)
        .join(GENOMAD.out.virus_summary)
        .join(SAMTOOLS_SORT.out.coverage)
        .join(BLASTN_VIROID.out.viroid)
    EVIDENCE(
        evidence_inputs,
        Channel.value(rvdb_taxmap),
        Channel.value(uniref90_taxmap),
        Channel.value(viroid_seqids),
        Channel.value(aggregate_script)
    )
    MULTIQC(
        FASTP.out.json.mix(SORTMERNA.out.log).collect()
    )
}
