// example: sbatch nextflow_wrapper.sbatch --data ont_subset --algorithm isoquant --result_name ont_subset/isoquant/test1

// input options

params.data = "isoseq" //options: --data "ont_subset", "ont", "isoseq", "isoseq_trimmed"
params.algorithm = "isoseq" //options: --algorithm "flair", "talon", "isoquant", "isoseq", "mandalorion"
params.use_sr = false
params.sr_config = "star" //options: --sr_config "jfs", "star", "star_p4", "star_p3", "star_p2", "star_old", "star_old_p3", "star_old_p2"
params.stringent = false
params.result_name = ""
params.joblog = ""
params.sqanti_filter = false
params.force_id_ignore = true

// paths
params.ont_metadata_samples_subset="/storage/gge/Fabian/nih/data/metadata/ont_samples_subset.tsv"
params.ont_metadata_concat_subset="/storage/gge/Fabian/nih/data/metadata/ont_concat_samples_subset.tsv"
params.ont_metadata_samples="/storage/gge/Fabian/nih/data/metadata/ont_samples.tsv"
params.ont_metadata_concat="/storage/gge/Fabian/nih/data/metadata/ont_concat_samples.tsv"
params.isoseq_metadata_samples="/storage/gge/Fabian/nih/data/metadata/isoseq_fl_samples.tsv"
params.isoseq_metadata_concat="/storage/gge/Fabian/nih/data/metadata/isoseq_fl_concat_samples.tsv"
params.isoseq_trimmed_metadata_samples="/storage/gge/Fabian/nih/data/metadata/isoseq_trimmed_samples.tsv"
params.isoseq_trimmed_metadata_concat="/storage/gge/Fabian/nih/data/metadata/isoseq_trimmed_concat_samples.tsv"
params.isoseq_trimmed_metadata_samples_subset="/storage/gge/Fabian/nih/data/metadata/isoseq_trimmed_samples_subset.tsv"
params.isoseq_trimmed_metadata_concat_subset="/storage/gge/Fabian/nih/data/metadata/isoseq_trimmed_concat_samples_subset.tsv"
params.isoseq_pipeline_metadata_concat="/storage/gge/Alejandro/nih/nextflow/isoseq_pipeline_concat_samples.tsv"
params.isoseq_pipeline_metadata_samples="/storage/gge/Alejandro/nih/nextflow/isoseq_pipeline_samples.tsv"
params.isoseq_pipeline_metadata_concat_subset="/storage/gge/Alejandro/nih/nextflow/subset_isoseq_pipeline_concat_samples.tsv"
params.isoseq_pipeline_metadata_samples_subset="/storage/gge/Alejandro/nih/nextflow/subset_isoseq_pipeline_samples.tsv"

// old mice
params.ont_old_metadata_samples="/storage/gge/Fabian/nih/data/metadata/old_mice/ont_old_samples.tsv"
params.ont_old_metadata_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/ont_old_concat_samples.tsv"
params.isoseq_old_metadata_samples="/storage/gge/Fabian/nih/data/metadata/old_mice/isoseq_old_fl_samples.tsv"
params.isoseq_old_metadata_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/isoseq_old_fl_concat_samples.tsv"
params.isoseq_trimmed_old_metadata_samples="/storage/gge/Fabian/nih/data/metadata/old_mice/isoseq_trimmed_old_samples.tsv"
params.isoseq_trimmed_old_metadata_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/isoseq_trimmed_old_concat_samples.tsv"

// partial joins
// ONT
params.ont_metadata_4samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/ont/ont_4samples.tsv"
params.ont_metadata_4samples_concat="Fabian/nih/data/metadata/partial_joins/ont/ont_concat_4samples.tsv"
params.ont_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/ont/ont_3samples.tsv"
params.ont_metadata_3samples_concat="Fabian/nih/data/metadata/partial_joins/ont/ont_concat_3samples.tsv"
params.ont_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/ont/ont_2samples.tsv"
params.ont_metadata_2samples_concat="Fabian/nih/data/metadata/partial_joins/ont/ont_concat_2samples.tsv"
// IsoSeq
params.isoseq_metadata_4samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_4samples.tsv"
params.isoseq_metadata_4samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_concat_4samples.tsv"
params.isoseq_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_3samples.tsv"
params.isoseq_metadata_3samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_concat_3samples.tsv"
params.isoseq_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_2samples.tsv"
params.isoseq_metadata_2samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq/isoseq_fl_concat_2samples.tsv"
// IsoSeq_trimmed
params.isoseq_trimmed_metadata_4samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_4samples.tsv"
params.isoseq_trimmed_metadata_4samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_concat_4samples.tsv"
params.isoseq_trimmed_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_3samples.tsv"
params.isoseq_trimmed_metadata_3samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_concat_3samples.tsv"
params.isoseq_trimmed_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_2samples.tsv"
params.isoseq_trimmed_metadata_2samples_concat="/storage/gge/Fabian/nih/data/metadata/partial_joins/isoseq_trimmed/isoseq_trimmed_concat_2samples.tsv"

// old mice partial joins
params.ont_old_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/ont/ont_old_2samples.tsv"
params.ont_old_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/ont/ont_old_3samples.tsv"
params.ont_old_metadata_2samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/ont/ont_old_concat_2samples.tsv"
params.ont_old_metadata_3samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/ont/ont_old_concat_3samples.tsv"
params.isoseq_old_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq/isoseq_old_fl_2samples.tsv"
params.isoseq_old_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq/isoseq_old_fl_3samples.tsv"
params.isoseq_old_metadata_2samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq/isoseq_old_fl_concat_2samples.tsv"
params.isoseq_old_metadata_3samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq/isoseq_old_fl_concat_3samples.tsv"
params.isoseq_trimmed_old_metadata_2samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq_trimmed/isoseq_trimmed_old_2samples.tsv"
params.isoseq_trimmed_old_metadata_3samples="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq_trimmed/isoseq_trimmed_old_3samples.tsv"
params.isoseq_trimmed_old_metadata_2samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq_trimmed/isoseq_trimmed_old_concat_2samples.tsv"
params.isoseq_trimmed_old_metadata_3samples_concat="/storage/gge/Fabian/nih/data/metadata/old_mice/partial_joins/isoseq_trimmed/isoseq_trimmed_old_concat_3samples.tsv"

// short reads via STAR output SJ.out.tab
params.sr_junctions_star="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_config.tsv"
params.sr_junctions_star_concat="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_concat_config.tsv"

// short reads via STAR output SJ.out.tab, partial joins
params.sr_junctions_star_p4="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_config_p4.tsv"
params.sr_junctions_star_concat_p4="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_concat_config_p4.tsv"
params.sr_junctions_star_p3="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_config_p3.tsv"
params.sr_junctions_star_concat_p3="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_concat_config_p3.tsv"
params.sr_junctions_star_p2="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_config_p2.tsv"
params.sr_junctions_star_concat_p2="/storage/gge/Fabian/nih/data/metadata/flair/star/flair_sr_junc_concat_config_p2.tsv"

// short reads via STAR output SJ.out.tab, old mice
params.sr_old_junctions_star="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_config.tsv"
params.sr_old_junctions_star_concat="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_concat_config.tsv"
params.sr_old_junctions_star_p3="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_config_p3.tsv"
params.sr_old_junctions_star_concat_p3="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_concat_config_p3.tsv"
params.sr_old_junctions_star_p2="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_config_p2.tsv"
params.sr_old_junctions_star_concat_p2="/storage/gge/Fabian/nih/data/metadata/flair/star_old/flair_sr_old_junc_concat_config_p2.tsv"

// short reads via FLAIR utility junctions_from_sam extracted junctions (not working)
params.sr_junctions_jfs="/storage/gge/Fabian/nih/data/metadata/flair/jfs/flair_sr_junctions.tsv"
params.sr_junctions_jfs_concat="/storage/gge/Fabian/nih/data/metadata/flair/jfs/flair_sr_junctions_concat.tsv"

// params.metadata_samples="/storage/gge/Fabian/nih/data/metadata/ont_samples_subset.tsv"
// params.metadata_concat="/storage/gge/Fabian/nih/data/metadata/ont_concat_samples_subset.tsv"

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
params.polyA = "/home/apadepe/polyA_site.bed"

// tools
params.base_tools_location = "$HOME/tools"
// params.base_tools_location = "$HOME/lr_pipelines"
params.tama_location = "${params.base_tools_location}/tama"
// params.sqanti_location = "${params.base_tools_location}/SQANTI3"
params.sqanti_location = "${params.base_tools_location}/SQANTI3_dev"
params.mandalorion_location = "${params.base_tools_location}/Mandalorion"

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
    // choose data type
    if ( params.data == "ont_subset" ) {
        input_ind = params.ont_metadata_samples_subset
        input_concat = params.ont_metadata_concat_subset

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont" ) {
        input_ind = params.ont_metadata_samples
        input_concat = params.ont_metadata_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont_p4" ) {
        input_ind = params.ont_metadata_4samples
        input_concat = params.ont_metadata_4samples_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont_p3" ) {
        input_ind = params.ont_metadata_3samples
        input_concat = params.ont_metadata_3samples_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont_p2" ) {
        input_ind = params.ont_metadata_2samples
        input_concat = params.ont_metadata_2samples_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "isoseq" ) {
        input_ind = params.isoseq_metadata_samples
        input_concat = params.isoseq_metadata_concat
        if ( params.algorithm == "isoseq" ) {
            input_ind = params.isoseq_pipeline_metadata_samples
            input_concat = params.isoseq_pipeline_metadata_concat
        }

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed" ) {
        input_ind = params.isoseq_trimmed_metadata_samples
        input_concat = params.isoseq_trimmed_metadata_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_subset" ) {
        input_ind = params.isoseq_trimmed_metadata_samples_subset
        input_concat = params.isoseq_trimmed_metadata_concat_subset

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_subset" ) {
        // input_ind = params.isoseq_metadata_samples_subset
        // input_concat = params.isoseq_metadata_concat_subset
        if ( params.algorithm == "isoseq" ) {
            input_ind = params.isoseq_pipeline_metadata_samples_subset
            input_concat = params.isoseq_pipeline_metadata_concat_subset
        }

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_p4" ) {
        input_ind = params.isoseq_metadata_4samples
        input_concat = params.isoseq_metadata_4samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_p3" ) {
        input_ind = params.isoseq_metadata_3samples
        input_concat = params.isoseq_metadata_3samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_p2" ) {
        input_ind = params.isoseq_metadata_2samples
        input_concat = params.isoseq_metadata_2samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_p4" ) {
        input_ind = params.isoseq_trimmed_metadata_4samples
        input_concat = params.isoseq_trimmed_metadata_4samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_p3" ) {
        input_ind = params.isoseq_trimmed_metadata_3samples
        input_concat = params.isoseq_trimmed_metadata_3samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_p2" ) {
        input_ind = params.isoseq_trimmed_metadata_2samples
        input_concat = params.isoseq_trimmed_metadata_2samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "ont_old" ) {
        input_ind = params.ont_old_metadata_samples
        input_concat = params.ont_old_metadata_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont_old_p3" ) {
        input_ind = params.ont_old_metadata_3samples
        input_concat = params.ont_old_metadata_3samples_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "ont_old_p2" ) {
        input_ind = params.ont_old_metadata_2samples
        input_concat = params.ont_old_metadata_2samples_concat

        data_type = "nanopore"
        fl_data = false
    }
    else if ( params.data == "isoseq_old" ) {
        input_ind = params.isoseq_old_metadata_samples
        input_concat = params.isoseq_old_metadata_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_old_p3" ) {
        input_ind = params.isoseq_old_metadata_3samples
        input_concat = params.isoseq_old_metadata_3samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_old_p2" ) {
        input_ind = params.isoseq_old_metadata_2samples
        input_concat = params.isoseq_old_metadata_2samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_old" ) {
        input_ind = params.isoseq_trimmed_old_metadata_samples
        input_concat = params.isoseq_trimmed_old_metadata_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_old_p3" ) {
        input_ind = params.isoseq_trimmed_old_metadata_3samples
        input_concat = params.isoseq_trimmed_old_metadata_3samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else if ( params.data == "isoseq_trimmed_old_p2" ) {
        input_ind = params.isoseq_trimmed_old_metadata_2samples
        input_concat = params.isoseq_trimmed_old_metadata_2samples_concat

        data_type = "pacbio_ccs"
        fl_data = true
    }
    else {
        error "Invalid data argument: ${params.data}"
    }

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

    // ** isoform characterization algorithm **
    // performs isoform characterization for each sample individually AND for concatenated sample (Join & Call) per tissue
    if ( params.algorithm == "flair" ) {

        if ( params.use_sr ) {
            if ( params.sr_config == "jfs" ) {
                sr_junctions=params.sr_junctions_jfs
                sr_junctions_concat=params.sr_junctions_jfs_concat
            }
            else if ( params.sr_config == "star" ) {
                sr_junctions=params.sr_junctions_star
                sr_junctions_concat=params.sr_junctions_star_concat
            }
            else if ( params.sr_config == "star_p4" ) {
                sr_junctions=params.sr_junctions_star_p4
                sr_junctions_concat=params.sr_junctions_star_concat_p4
            }
            else if ( params.sr_config == "star_p3" ) {
                sr_junctions=params.sr_junctions_star_p3
                sr_junctions_concat=params.sr_junctions_star_concat_p3
            }
            else if ( params.sr_config == "star_p2" ) {
                sr_junctions=params.sr_junctions_star_p2
                sr_junctions_concat=params.sr_junctions_star_concat_p2
            }
            else if ( params.sr_config == "star_old" ) {
                sr_junctions=params.sr_old_junctions_star
                sr_junctions_concat=params.sr_old_junctions_star_concat
            }
            else if ( params.sr_config == "star_old_p3" ) {
                sr_junctions=params.sr_old_junctions_star_p3
                sr_junctions_concat=params.sr_old_junctions_star_concat_p3
            }
            else if ( params.sr_config == "star_old_p2" ) {
                sr_junctions=params.sr_old_junctions_star_p2
                sr_junctions_concat=params.sr_old_junctions_star_concat_p2
            } else {
                error "Invalid sr_config argument: ${params.sr_config}"
            }
        }
        else {
            sr_junctions=params.sr_junctions_star
            sr_junctions_concat=params.sr_junctions_star_concat
        }

        algo_output = flair(input_ind, sr_junctions, sr_junctions_concat, wd, joblog)
    }
    else if ( params.algorithm == "talon" ) {
        algo_output = talon(input_ind, input_concat, wd, joblog)
    }
    else if ( params.algorithm == "isoquant" ) {
        algo_output = isoquant(input_ind, input_concat, data_type, fl_data, wd, joblog)
    }
    else if ( params.algorithm == "mandalorion" ) {
        algo_output = mandalorion(input_ind, input_concat, data_type, wd, joblog)
    }
    else if ( params.algorithm == "isoseq" ) {
        if ( data_type == "pacbio_ccs") {
            params.force_id_ignore = false
            algo_output = isoseq(input_ind, input_concat, wd, joblog)
            if ( params.sqanti_filter ) {
                sqanti_orth_output = sqanti3_orthogonal(algo_output.metadata_ind, algo_output.metadata_concat, wd, joblog)
                algo_output = sqanti3_filter(sqanti_orth_output.metadata_orth_ind, sqanti_orth_output.metadata_orth_concat, wd, joblog)
            }
            algo_output = isoseq_metadata(algo_output.metadata_ind, algo_output.metadata_concat, wd)
        } else {
            error "Invalid data argument with isoseq pipeline: ${params.data}"
        }
    }
    else if ( params.algorithm == "bambu" ){
        algo_output = bambu(input_ind, input_concat, wd, joblog)
    }
    else {
        error "Invalid algorithm argument: ${params.algorithm}"
    }

    // ** Call & Join vs. Join & Call comparison **
    // run TAMA to combine results of individual samples
    tama_condition_output = tama_condition(algo_output.metadata_ind, wd, joblog)

    // run TAMA across all conditions to get uniform transcript ids
    tama_full_output = tama_full(algo_output.metadata_ind, algo_output.metadata_concat, wd, joblog)

    // combine the count matrices of the individual samples
    merge_counts_output = merge_counts(algo_output.metadata_ind, algo_output.metadata_concat, tama_condition_output.metadata_tama_condition, tama_full_output.metadata_tama_full, wd, joblog)

    // run SQANTI3 quality control on all gtfs
    sq3_output = sqanti3(merge_counts_output.metadata_merged, wd, joblog)

    // prepare data for report
    prepare_report(report, sq3_output.metadata_sq3, wd)
    // after this, download the folder and locally run the report to create the plots
}
