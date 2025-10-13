// example: sbatch nextflow_wrapper.sbatch --data ont_subset --algorithm isoquant --result_name ont_subset/isoquant/test1

// input options

params.data = "isoseq" //options: --data "ont_subset", "ont", "isoseq", "isoseq_trimmed"
params.algorithm = "isoseq" //options: --algorithm "flair", "talon", "isoquant", "isoseq", "mandalorion"
params.use_sr = false
params.sr_config = "star" //options: --sr_config "jfs", "star", "star_p4", "star_p3", "star_p2", "star_old", "star_old_p3", "star_old_p2"
params.stringent = false
params.result_name = ""
params.joblog = ""
params.force_id_ignore = true
params.algo_metadata_ind = ""
params.algo_metadata_concat = ""

params.data_dir = "${projectDir}/../../../data/output"
params.report_dir = "${projectDir}/../../../reports"
params.empty_report = "/storage/gge/Fabian/nih/empty_report"
params.src_dir = "${projectDir}/../scripts"
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
params.tama_location = "${params.base_tools_location}/tama"
// params.sqanti_location = "${params.base_tools_location}/SQANTI3"
params.sqanti_location = "${params.base_tools_location}/SQANTI3_dev"
params.mandalorion_location = "${params.base_tools_location}/Mandalorion"
params.flair_location = "${params.base_tools_location}/flair"

// include statements need to be AFTER parameter definitions
include { flair } from '../modules/algorithms/flair.nf'
include { talon } from '../modules/algorithms/talon.nf'
include { isoquant } from '../modules/algorithms/isoquant.nf'
include { bambu } from '../modules/algorithms/bambu.nf'
include { mandalorion } from '../modules/algorithms/mandalorion.nf'
include { isoseq; isoseq_metadata } from '../modules/algorithms/isoseq.nf'
include { tama_condition; tama_full } from '../modules/tama.nf'
include { sqanti3; sqanti3_orthogonal; sqanti3_filter } from '../modules/sqanti3.nf'
include { merge_counts } from '../modules/merge_counts.nf'
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
        joblog = "${report}/joblog.tsv"
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


    // ** Call & Join vs. Join & Call comparison **
    // run TAMA to combine results of individual samples
    tama_condition_output = tama_condition(params.algo_metadata_ind, wd, joblog)

    // run TAMA across all conditions to get uniform transcript ids
    tama_full_output = tama_full(params.algo_metadata_ind, params.algo_metadata_concat, wd, joblog)

    // combine the count matrices of the individual samples
    merge_counts_output = merge_counts(params.algo_metadata_ind, params.algo_metadata_concat, tama_condition_output.metadata_tama_condition, tama_full_output.metadata_tama_full, wd, joblog)

    // run SQANTI3 quality control on all gtfs
    sq3_output = sqanti3(merge_counts_output.metadata_merged, wd, joblog)

    // prepare data for report
    prepare_report(report, sq3_output.metadata_sq3, wd)
    // after this, download the folder and locally run the report to create the plots
}
