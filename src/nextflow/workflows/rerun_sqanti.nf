// example: sbatch nf_wrapper_rerun_sqanti.sbatch --metadata_merged /storage/gge/Alejandro/for_fabian/paper_stuff/isoseq_bambu_run1_data/merge_counts/metadata_merged.tsv --result_name alejandro/isoseq/bambu/run1

// input options

params.result_name = ""
params.joblog = ""
params.sqanti_filter = false
params.force_id_ignore = true

params.data_dir = "${projectDir}/../../../data/output"
params.report_dir = "${projectDir}/../../../reports"
params.empty_report = "/storage/gge/Fabian/nih/empty_report"
params.src_dir = "${projectDir}/../scripts"
params.metadata_merged = ""
// params.report_dir = "/home/apadepe/documenting_NIH/fabian/reports"
// params.empty_report = "/storage/gge/Fabian/nih/empty_report"
// params.src_dir = "/home/apadepe/documenting_NIH/fabian/src/nextflow/scripts"

// files
params.genome = "/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
params.annotation = "/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
params.annotation_db = "/storage/gge/Fabian/nih/data/metadata/isoquant/mm39.ncbiRefSeq_SIRV.db"
params.cage = "/storage/gge/genomes/mouse_ref_NIH/reference_genome/lft_mm39_CAGE.bed"

// tools
params.base_tools_location = "$HOME/tools"
// params.base_tools_location = "$HOME/lr_pipelines"
// params.sqanti_location = "${params.base_tools_location}/SQANTI3"
params.sqanti_location = "${params.base_tools_location}/SQANTI3_dev"

// include statements need to be AFTER parameter definitions
include { sqanti3; sqanti3_orthogonal; sqanti3_filter } from '../modules/sqanti3.nf'
include { prepare_report } from '../modules/prepare_report.nf'

workflow {

    if (params.result_name == null || params.result_name.isEmpty()) {
        wd = "${params.data_dir}/${params.data}/${params.algorithm}_data"
        report = "${params.report_dir}/${params.data}/${params.algorithm}_report"
    } else {
        wd = "${params.data_dir}/${params.result_name}_data"
        report = "${params.report_dir}/${params.result_name}_report"
    }

    if (params.joblog == null || params.joblog.isEmpty()) {
        joblog = "${report}/joblog_sq.tsv"
    } else {
        joblog = "${report}/${params.joblog}"
    }

    println "Working directory: ${wd}"
    println "Report directory: ${report}"
    println "Job log: ${joblog}"

    // Check if the working directories exist, if not, create them
    dataDir = new File(wd)
    if (!dataDir.exists()) {
        dataDir.mkdirs()
    }
    reportDir = new File(report)
    if (!reportDir.exists()) {
        reportDir.mkdirs()
    }

    // record the job id of the workflow process
    new File(joblog).withWriter { writer ->
        writer.writeLine("WORKFLOW\t${System.getenv('SLURM_JOB_ID')}")
    }

    // run SQANTI3 quality control on all gtfs
    sq3_output = sqanti3(params.metadata_merged, wd, joblog)

    // prepare data for report
    prepare_report(report, sq3_output.metadata_sq3, wd)
    // after this, download the folder and locally run the report to create the plots
}
