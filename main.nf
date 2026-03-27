nextflow.enable.dsl=2

// --- Import Modules ---
include { SHEETMAKER } from './modules/preprocessing'
include { TRIMGALORE } from './modules/trimgalore'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIM; FASTQC as FASTQC_RIBO } from './modules/fastqc'
include { MULTIQC as MULTIQC_RAW; MULTIQC as MULTIQC_TRIM; MULTIQC as MULTIQC_RIBO } from './modules/multiqc'
include { RIBODETECTOR } from './modules/ribodetector_Singularity'
include { STAR_ALIGN; FEATURECOUNTS; MERGE_COUNTS } from './modules/alignment'
include { KRAKEN2 } from './modules/kraken2'
include { BRACKEN } from './modules/bracken'
include { KRONA } from './modules/krona'
include { SALMON_INDEX; SALMON_QUANT } from './modules/salmon'
include { DGE_R_ANALYSIS } from './modules/dge'

// Added (were previously never called)
include { SORTMERNA } from './modules/sortmerna'
include { PARABRICKS_RNA } from './modules/parabricks'

// --- DGE-only Python module (NEW) ---
include { DGE_ANALYSIS } from './modules/dge_analysis'

// --- Main Workflow ---
workflow {

    // If the user provides precomputed counts + metadata, run DGE without requiring reads.
    if (params.run_dge && params.dge_counts && params.dge_metadata) {
        DGE_R_ANALYSIS(
            file(params.dge_counts),
            file(params.dge_metadata),
            params.dge_control,
            params.dge_treatment,
            params.dge_comparison_name,
            params.dge_tool,
            params.batch_method,
            params.batch_col,
            params.group_col,
            params.covariates
        )
        return
    }

    // 1. Data Discovery
    if (!params.reads) {
        error "Please provide --reads parameter (or run DGE-only by providing --run_dge true --dge_counts <file> --dge_metadata <file>)"
    }
    def input_path = file(params.reads).toString()
    samplesheet_ch = SHEETMAKER(input_path)

    samplesheet_ch.csv
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, [file(row.fastq_1), file(row.fastq_2)]) }
        .set { read_pairs }

    // 2. Initial QC and Trimming
    if (params.run_trimgalore) {
        if (!params.skip_fastqc_raw) {
            FASTQC_RAW(read_pairs, "_raw")
            if (!params.skip_multiqc_raw) {
                MULTIQC_RAW(FASTQC_RAW.out.zip.collect(), "_raw")
            }
        }

        // 3. Trimming
        TRIMGALORE(read_pairs)
        ch_after_trim = TRIMGALORE.out.reads

        if (!params.skip_fastqc_trimmed) {
            FASTQC_TRIM(TRIMGALORE.out.reads, "_trimmed")
            if (!params.skip_multiqc_trimmed) {
                MULTIQC_TRIM(FASTQC_TRIM.out.zip.collect(), "_trimmed")
            }
        }
    } else {
        ch_after_trim = read_pairs
    }

    // 3.6. Kraken2 & Bracken Taxonomy (optional) -- run after Trimming
    if (params.run_kraken2) {
        KRAKEN2(ch_after_trim, file(params.kraken2_db))

        // Run BRACKEN on Kraken2 reports
        BRACKEN(KRAKEN2.out.report, file(params.kraken2_db))

        // Run KRONA on Kraken2 reports
        KRONA(KRAKEN2.out.report)
    }

    // 4. RiboDetector (GPU) (optional) -- run BEFORE SortMeRNA
    ch_after_ribo = ch_after_trim
    if (params.run_ribodetector) {
        RIBODETECTOR(ch_after_trim)
        ch_after_ribo = RIBODETECTOR.out.reads

        if (!params.skip_fastqc_ribo) {
            FASTQC_RIBO(RIBODETECTOR.out.reads, "_ribo")
            if (!params.skip_multiqc_ribo) {
                MULTIQC_RIBO(FASTQC_RIBO.out.zip.collect(), "_ribo")
            }
        }
    }

    // 3.5. SortMeRNA (optional) -- run AFTER RiboDetector
    ch_to_use = ch_after_ribo
    if (params.run_sortmerna) {
        SORTMERNA(ch_after_ribo, file(params.smr_db_dir))
        ch_to_use = SORTMERNA.out.reads
    }

    // 5. Parabricks (optional, GPU)
    if (params.run_parabricks) {
        PARABRICKS_RNA(
            ch_to_use,
            file(params.parabricks_index),
            file(params.genome_fasta),
            file(params.genome_gtf)
        )
    }

    // 6. Alignment and Quantification
    ch_fc_counts = Channel.empty()
    ch_star_counts = Channel.empty()
    ch_salmon_counts = Channel.empty()

    if (params.run_star) {
        STAR_ALIGN(
            ch_to_use,
            file(params.star_index),
            file(params.genome_gtf)
        )
        ch_star_counts = STAR_ALIGN.out.counts.collect()

        if (params.run_featurecounts) {
            FEATURECOUNTS(
                STAR_ALIGN.out.bam,
                file(params.genome_gtf)
            )
            ch_fc_counts = FEATURECOUNTS.out.counts
                        .map { sample_id, fc_file -> fc_file }
                        .collect()
        }
    }

    if (params.run_salmon) {
        if (params.salmon_index) {
            ch_salmon_index = file(params.salmon_index)
        } else {
            SALMON_INDEX(file(params.genome_fasta))
            ch_salmon_index = SALMON_INDEX.out.index
        }
        SALMON_QUANT(
            ch_to_use,
            ch_salmon_index,
            file(params.genome_gtf)
        )
        ch_salmon_counts = SALMON_QUANT.out.quant.collect()
    }

    if (params.run_star || params.run_salmon || params.run_featurecounts) {
        MERGE_COUNTS(
            ch_fc_counts.ifEmpty([]),
            ch_star_counts.ifEmpty([]),
            ch_salmon_counts.ifEmpty([])
        )

        if (params.run_dge) {
            DGE_R_ANALYSIS(
                MERGE_COUNTS.out.matrix,
                file(params.samplesheet), // Using samplesheet as metadata or separate metadata param
                params.dge_control,
                params.dge_treatment,
                params.dge_comparison_name,
                params.dge_tool,
                params.batch_method,
                params.batch_col,
                params.group_col,
                params.covariates
            )
        }
    }
}
