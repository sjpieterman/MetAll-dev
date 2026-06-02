nextflow.enable.dsl=2

include { KRAKEN2 } from '../modules/kraken2'
include { BRACKEN } from '../modules/bracken'
include { KRONA } from '../modules/krona'
include { IDENTIFY_CONTAMINANTS; APPLY_DECONTAM } from '../modules/decontam'

workflow RUN_TAXONOMY {
    take:
    reads_ch
    suffix

    main:
    // Append suffix to sample ID for unique filenames
    ch_input = reads_ch.map { id, reads -> tuple("${id}${suffix}", reads) }
    
    KRAKEN2(ch_input, file(params.kraken2_db))
    
    ch_kraken_reports = KRAKEN2.out.report
    
    if (params.run_decontam && params.decontam_metadata) {
        IDENTIFY_CONTAMINANTS(
            ch_kraken_reports.map { it[1] }.collect(),
            file(params.decontam_metadata),
            params.decontam_col,
            params.decontam_method,
            params.decontam_threshold
        )
        APPLY_DECONTAM(ch_kraken_reports, IDENTIFY_CONTAMINANTS.out.contaminants)
        
        // Fallback: Use decontaminated report if available, otherwise original Kraken2 report
        ch_reports = ch_kraken_reports
            .join(APPLY_DECONTAM.out.report, remainder: true)
            .map { id, kraken, decontam ->
                [ id, decontam ?: kraken ]
            }
    } else {
        ch_reports = ch_kraken_reports
    }

    BRACKEN(ch_reports, file(params.kraken2_db))
    KRONA(ch_reports)
    
    emit:
    reports = ch_reports
}
