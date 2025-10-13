process isoquant {
    input:
    path in_metadata_samples
    path in_metadata_concat
    val data_type
    val fl_data
    path wd
    path joblog

    output:
    path "isoquant_out_metadata_ind.tsv", emit: metadata_ind
    path "isoquant_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    echo "data_type: ${data_type}"
    echo "fl_data: ${fl_data}"

    ${params.src_dir}/algorithms/isoquant/run_isoquant.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation_db}" \\
        --metadata_samples "${in_metadata_samples}" \\
        --metadata_concat "${in_metadata_concat}" \\
        --data_type "${data_type}" \\
        --fl_data ${fl_data} \\
        --joblog ${joblog}

    ln -s "\$PWD" "${wd}/isoquant"
    """
}
