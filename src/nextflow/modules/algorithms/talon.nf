process talon {
    input:
    path in_metadata_samples
    path in_metadata_concat
    path wd
    path joblog

    output:
    path "talon_out_metadata_ind.tsv", emit: metadata_ind
    path "talon_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/algorithms/talon/run_talon.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_samples "${in_metadata_samples}" \\
        --metadata_concat "${in_metadata_concat}" \\
        --joblog "${joblog}"${params.talon_reformat_skip_gene_id_gene_name ? ' --skip_gene_id_gene_name' : ''}${params.talon_max_reads > 0 ? " --talon_max_reads ${params.talon_max_reads}" : ''}

    ln -s "\$PWD" ${wd}/talon
    """
}
