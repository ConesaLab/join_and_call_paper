process mandalorion {
    input:
    path in_metadata_samples
    path in_metadata_concat
    val data_type
    path wd
    path joblog

    output:
    path "mandalorion_out_metadata_ind.tsv", emit: metadata_ind
    path "mandalorion_out_metadata_concat.tsv", emit: metadata_concat

    beforeScript 'mkdir log'

    script:
    """
    echo "data_type: ${data_type}"
    if ["${data_type}" != "pacbio_ccs"]; then
        echo "Mandalorion should only be used on PacBio or R2C2 reads"
        exit 1
    fi

    ${params.src_dir}/algorithms/mandalorion/run_mandalorion.sh \\
        --wd "." \\
        --genome "${params.genome}" \\
        --annotation "${params.annotation}" \\
        --metadata_samples "${in_metadata_samples}" \\
        --metadata_concat "${in_metadata_concat}" \\
        --mandalorion_location "${params.mandalorion_location}" \\
        --joblog ${joblog}

    ln -s "\$PWD" "${wd}/mandalorion"
    """
}
