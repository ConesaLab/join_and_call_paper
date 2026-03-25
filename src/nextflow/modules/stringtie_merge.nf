process stringtie_merge_condition {
    input:
    path metadata_ind
    path wd
    path joblog

    output:
    path "metadata_stmerge_condition.tsv", emit: metadata_stmerge_condition

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/stringtie_merge/run_stringtie_merge_condition.sh \\
        --wd "." \\
        --metadata_samples "${metadata_ind}" \\
        --joblog "${joblog}"

    ln -snf "\$PWD" ${wd}/stringtie_merge_condition
    """
}

process sqanti3_stmerge {
    input:
    path metadata_stmerge
    path wd
    path joblog

    output:
    path "metadata_sq3_stmerge.tsv", emit: metadata_sq3_stmerge

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/stringtie_merge/run_sqanti3_stmerge.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_stmerge "${metadata_stmerge}" \\
        --sqanti_path "${params.sqanti_location}" \\
        --joblog "${joblog}"

    ln -snf "\$PWD" ${wd}/sqanti3_stmerge
    """
}
