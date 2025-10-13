process isoseq {
    input:
    path metadata_samples
    path metadata_concat
    path wd
    path joblog


    output:
    path "isoseq_out_metadata_ind.tsv", emit: metadata_ind
    path "isoseq_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    ${params.src_dir}/algorithms/isoseq/run_isoseq.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --metadata_samples "${metadata_samples}" \\
        --metadata_concat "${metadata_concat}" \\
        --joblog "${joblog}"
        
    ln -s "\$PWD" ${wd}/isoseq
    """
}

process isoseq_metadata {
    input:
    path metadata_samples
    path metadata_concat
    path wd

    output:
    path "isoseq_final_metadata.tsv", emit: metadata_ind
    path "isoseq_final_metadata_concat.tsv", emit: metadata_concat

    script:
    """
    cut -f1,7,8 ${metadata_concat} > "isoseq_final_metadata_concat.tsv"
    cut -f1-5,7,8 ${metadata_samples} > "isoseq_final_metadata.tsv"

    ln -s "\$PWD" ${wd}/isoseq_metadata
    """
}