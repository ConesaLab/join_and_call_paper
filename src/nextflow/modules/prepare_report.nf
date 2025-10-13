process prepare_report{
    input:
    path target
    path metadata_merged
    path wd

    shell:
    '''
    !{params.src_dir}/prepare_report/prepare_report.sh \\
        --wd "!{target}" \\
        --empty_report "!{params.empty_report}" \\
        --genome "!{params.genome}" \\
        --metadata_merged "!{metadata_merged}"

    ln -s "\$PWD" !{wd}/prepare_report
    ln -s "\$(realpath !{target})" !{wd}/report
    '''
}