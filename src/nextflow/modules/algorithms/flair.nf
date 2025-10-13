process flair {
    input:
    path metadata_samples
    path sr_junctions
    path sr_junctions_concat
    path wd
    path joblog

    output:
    path "flair_out_metadata_ind.tsv", emit: metadata_ind
    path "flair_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/algorithms/flair/run_flair.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_samples "${metadata_samples}" \\
        --use_sr "${params.use_sr}" \\
        --stringent "${params.stringent}" \\
        --sr_junctions "${sr_junctions}" \\
        --sr_junctions_concat "${sr_junctions_concat}" \\
        --joblog "${joblog}"

    ln -s "\$PWD" ${wd}/flair
    """
}
