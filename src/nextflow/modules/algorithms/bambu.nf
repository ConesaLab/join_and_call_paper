process bambu {
    input:
    path in_metadata_samples
    path in_metadata_concat
    path wd
    path joblog

    output:
    path "bambu_out_metadata_ind.tsv", emit: metadata_ind
    path "bambu_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/algorithms/bambu/run_bambu.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_samples "${in_metadata_samples}" \\
        --metadata_concat "${in_metadata_concat}" \\
        --joblog "${joblog}" \\

    ln -s "\$PWD" ${wd}/bambu
    """
}
