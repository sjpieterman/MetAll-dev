nextflow.enable.dsl=2

// --- Import Modules ---
include { SHEETMAKER } from './modules/preprocessing'
include { TRIMGALORE } from './modules/trimgalore'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIM; FASTQC as FASTQC_RIBO; FASTQC as FASTQC_RIBO_RRNA; FASTQC as FASTQC_SMR; FASTQC as FASTQC_SMR_RRNA } from './modules/fastqc'
include { MULTIQC as MULTIQC_RAW; MULTIQC as MULTIQC_TRIM; MULTIQC as MULTIQC_RIBO; MULTIQC as MULTIQC_SMR } from './modules/multiqc'
include { RIBODETECTOR } from './modules/ribodetector_Singularity'
include { STAR_ALIGN; FEATURECOUNTS; MERGE_COUNTS } from './modules/alignment'
include { SALMON_INDEX; SALMON_QUANT } from './modules/salmon'
include { SALMON_VIRULENCE_INDEX; SALMON_VIRULENCE_QUANT } from './modules/virulence'
include { DGE_R_ANALYSIS } from './modules/dge'
include { FASTP } from './modules/fastp'

// Added (were previously never called)
include { SORTMERNA } from './modules/sortmerna'
include { PARABRICKS_RNA } from './modules/parabricks'

// --- DGE-only Python module (NEW) ---
include { DGE_ANALYSIS } from './modules/dge_analysis'

include { RUN_TAXONOMY as RUN_TAXONOMY_RRNA; RUN_TAXONOMY as RUN_TAXONOMY_UNMAPPED } from './subworkflows/taxonomy'

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
    if (params.run_qc) {
        FASTQC_RAW(read_pairs, "_raw")
        MULTIQC_RAW(FASTQC_RAW.out.zip.collect(), "_raw")
    }

    if (params.use_fastp) {
        // Fastp handles both QC and Trimming in one fast pass
        FASTP(read_pairs)
        ch_after_trim = FASTP.out.reads

        if (params.run_qc) {
            FASTQC_TRIM(FASTP.out.reads, "_trimmed")
            MULTIQC_TRIM(FASTP.out.json.mix(FASTQC_TRIM.out.zip).collect(), "_trimmed")
        }
    } else if (params.run_trimgalore) {
        // 3. Trimming
        TRIMGALORE(read_pairs)
        ch_after_trim = TRIMGALORE.out.reads

        if (params.run_qc) {
            FASTQC_TRIM(TRIMGALORE.out.reads, "_trimmed")
            MULTIQC_TRIM(FASTQC_TRIM.out.zip.collect(), "_trimmed")
        }
    } else {
        ch_after_trim = read_pairs
    }


    // 4. RiboDetector (GPU) (optional) -- run BEFORE SortMeRNA
    ch_after_ribo = ch_after_trim
    ch_rrna_for_taxonomy = Channel.empty()
    if (params.run_ribodetector) {
        RIBODETECTOR(ch_after_trim)
        ch_after_ribo = RIBODETECTOR.out.reads
        ch_rrna_for_taxonomy = RIBODETECTOR.out.rrna_reads

        if (params.run_qc) {
            FASTQC_RIBO(RIBODETECTOR.out.reads, "_ribo_free")
            FASTQC_RIBO_RRNA(RIBODETECTOR.out.rrna_reads, "_ribo_rrna")
            MULTIQC_RIBO(
                FASTQC_RIBO.out.zip
                .mix(FASTQC_RIBO_RRNA.out.zip)
                .collect(), 
                "_ribo"
            )
        }
    }

    // 3.5. SortMeRNA (optional) -- run AFTER RiboDetector
    ch_to_use = ch_after_ribo
    if (params.run_sortmerna) {
        SORTMERNA(ch_after_ribo, file(params.smr_db_dir))
        ch_to_use = SORTMERNA.out.reads
        ch_rrna_for_taxonomy = SORTMERNA.out.rrna_reads

        if (params.run_qc) {
            FASTQC_SMR(SORTMERNA.out.reads, "_smr_free")
            FASTQC_SMR_RRNA(SORTMERNA.out.rrna_reads, "_smr_rrna")
            MULTIQC_SMR(
                FASTQC_SMR.out.zip
                .mix(FASTQC_SMR_RRNA.out.zip)
                .collect(), 
                "_smr"
            )
        }
    }

    // 4.5 Taxonomy on rRNA (optional)
    if (params.run_kraken2) {
        RUN_TAXONOMY_RRNA(ch_rrna_for_taxonomy, ".rRNA")
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

    // 5.5 Virulence Factor Analysis (optional)
    if (params.run_virulence) {
        if (params.virulence_index && params.virulence_index != true && params.virulence_index != "null") {
            ch_virulence_index = file(params.virulence_index)
        } else if (params.virulence_db && params.virulence_db != true && params.virulence_db != "null") {
            SALMON_VIRULENCE_INDEX(file(params.virulence_db))
            ch_virulence_index = SALMON_VIRULENCE_INDEX.out.index
        } else {
            error "Please provide --virulence_db or --virulence_index for virulence factor analysis"
        }
        SALMON_VIRULENCE_QUANT(
            ch_to_use,
            ch_virulence_index
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

        if (params.run_kraken2) {
            RUN_TAXONOMY_UNMAPPED(STAR_ALIGN.out.unmapped, ".unmapped")
        }

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
        if (params.salmon_index && params.salmon_index != true && params.salmon_index != "null") {
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
                params.samplesheet ? file(params.samplesheet) : samplesheet_ch.csv,
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
