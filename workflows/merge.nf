
process merge_calls {
    label "wf_basecalling"
    cpus 4
    input:
        path(ref)
        path(crams)
        val(filetag)
    output:
        tuple path("${params.sample_name}.${filetag}.cram"), path("${params.sample_name}.${filetag}.cram.crai")
    script:
    """
    samtools merge ${params.sample_name}.${filetag}.cram ${crams} --no-PG -O CRAM --write-index --reference ${ref} --threads ${task.cpus}
    """
}
