process sqanti3{
    input:
    path metadata_merged
    path wd
    path joblog

    output:
    path "metadata_sq3.tsv", emit: metadata_sq3

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/sqanti3/run_sqanti3_qc.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_merged "${metadata_merged}" \\
        --sqanti_path "${params.sqanti_location}" \\
        --joblog "${joblog}"

    ln -s "\$PWD" ${wd}/sqanti3
    """
}

process sqanti3_orthogonal {
    input:
    path metadata_ind
    path metadata_concat
    path wd
    path joblog

    output:
    path "metadata_orth_ind.tsv", emit: metadata_orth_ind
    path "metadata_orth_concat.tsv", emit: metadata_orth_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/sqanti3/run_sqanti3_qc_orthogonal.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_ind "${metadata_ind}" \\
        --metadata_concat "${metadata_concat}" \\
        --cage "${params.cage}" \\
        --polyA "${params.polyA}" \\
        --sqanti_path "${params.sqanti_location}" \\
        --joblog "${joblog}"
        
    ln -s "\$PWD" ${wd}/sqanti_orth
    """
}

process sqanti3_filter {
    input:
    path metadata_orth_ind
    path metadata_orth_concat
    path wd
    path joblog

    output:
    path "metadata_filter_ind.tsv", emit: metadata_ind
    path "metadata_filter_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/sqanti3/run_sqanti3_filter.sh \\
        --wd "." \\
        --metadata_ind "${metadata_orth_ind}" \\
        --metadata_concat "${metadata_orth_concat}" \\
        --sqanti_path "${params.sqanti_location}" \\
        --joblog "${joblog}"
        
    ln -s "\$PWD" ${wd}/sqanti_filter
    """
}