process tama_condition {
    input:
    path metadata_ind
    path wd
    path joblog

    output:
    path "metadata_tama_condition.tsv", emit: metadata_tama_condition

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/tama/run_tama_condition.sh \\
        --wd "." \\
        --metadata_samples "${metadata_ind}" \\
        --tama_path "${params.tama_location}" \\
        --joblog "${joblog}"

    ln -s "\$PWD" ${wd}/tama_condition
    """
}

process tama_full {
    input:
    path metadata_ind
    path metadata_concat
    path wd
    path joblog

    output:
    path "metadata_tama_full.tsv", emit: metadata_tama_full

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/tama/run_tama_full.sh \\
        --wd "." \\
        --metadata_samples "${metadata_ind}" \\
        --metadata_concat "${metadata_concat}" \\
        --tama_path "${params.tama_location}" \\
        --joblog "${joblog}"

    ln -s "\$PWD" ${wd}/tama_full
    """
}