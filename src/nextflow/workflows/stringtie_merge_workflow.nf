// Standalone workflow: StringTie merge (per condition) + SQANTI3 QC
// Operates on existing pipeline results (read-only) and writes to a new output directory.
//
// Usage:
//   sbatch nf_wrapper_stringtie_merge.sbatch \
//       --metadata_ind /path/to/existing/metadata_ind.tsv \
//       --result_name isoseq/isoquant/stringtie_merge

params.metadata_ind = ""
params.result_name = ""
params.joblog = ""

params.data_dir = "${projectDir}/../../../data/output"
params.report_dir = "${projectDir}/../../../reports"
params.src_dir = "${projectDir}/../scripts"

params.genome = "/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
params.annotation = "/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"

params.base_tools_location = "$HOME/tools"
params.sqanti_location = "${params.base_tools_location}/SQANTI3_dev"

include { stringtie_merge_condition; sqanti3_stmerge } from '../modules/stringtie_merge.nf'

workflow {
    if (params.metadata_ind == null || params.metadata_ind.isEmpty()) {
        error "Please provide --metadata_ind (path to an existing algorithm metadata_ind.tsv)"
    }
    if (params.result_name == null || params.result_name.isEmpty()) {
        error "Please provide --result_name (e.g. isoseq/isoquant/stringtie_merge)"
    }

    wd = "${params.data_dir}/${params.result_name}_data"
    report = "${params.report_dir}/${params.result_name}_report"

    if (params.joblog == null || params.joblog.isEmpty()) {
        joblog = "${report}/joblog_stmerge.tsv"
    } else {
        joblog = "${report}/${params.joblog}"
    }

    println "Working directory: ${wd}"
    println "Report directory: ${report}"
    println "Job log: ${joblog}"

    dataDir = new File(wd)
    if (!dataDir.exists()) {
        dataDir.mkdirs()
    }
    reportDir = new File(report)
    if (!reportDir.exists()) {
        reportDir.mkdirs()
    }

    new File(joblog).withWriter { writer ->
        writer.writeLine("WORKFLOW\t${System.getenv('SLURM_JOB_ID')}")
    }

    stmerge_output = stringtie_merge_condition(params.metadata_ind, wd, joblog)

    sqanti3_stmerge(stmerge_output.metadata_stmerge_condition, wd, joblog)
}
