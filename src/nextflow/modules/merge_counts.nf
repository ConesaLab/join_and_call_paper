process merge_counts{
    input:
    path metadata_ind
    path metadata_concat
    path metadata_tama
    path metadata_tama_full
    path wd
    path joblog

    output:
    path "metadata_merged.tsv", emit: metadata_merged

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/merge_counts/merge_counts.sh \\
        --wd "." \\
        --metadata_ind "${metadata_ind}" \\
        --metadata_concat "${metadata_concat}" \\
        --metadata_tama "${metadata_tama}" \\
        --metadata_tama_full "${metadata_tama_full}" \\
        --joblog "${joblog}"

    ln -s "\$PWD" ${wd}/merge_counts
    """
}