
// this process is shared by both the uCRAM and CRAM arms of the basecalling workflow
// for uCRAM the staged ref is the OPTIONAL_FILE, so we withhold the ref arg
process merge_calls {
    label "wf_basecalling"
    cpus params.merge_threads
    input:
        path(ref)
        path(crams)
        val(filetag)
        tuple val(align_ext), val(index_ext) // either [bam, bai] or [cram, crai]
    output:
        tuple path("${params.sample_name}.${filetag}.${align_ext}"), path("${params.sample_name}.${filetag}.${align_ext}.${index_ext}")
    script:
    def ref_arg = ref.name != "OPTIONAL_FILE" ? "--reference ${ref}" : ""
    """
    samtools merge "${params.sample_name}.${filetag}.${align_ext}##idx##${params.sample_name}.${filetag}.${align_ext}.${index_ext}" ${crams} --no-PG -O ${align_ext} --write-index ${ref_arg} --threads ${task.cpus}
    """
}

process merge_calls_to_fastq {
    label "wf_basecalling"
    cpus { params.merge_threads + params.ubam_bam2fq_threads }
    input:
        path(crams)
        val(filetag)
    output:
        path("${params.sample_name}.${filetag}.fq.gz")
    script:
    """
    samtools merge ${crams} --no-PG -O CRAM -@ ${params.merge_threads} -o - | samtools bam2fq -T 1 -@ ${params.ubam_bam2fq_threads} -0 ${params.sample_name}.${filetag}.fq.gz -
    """
}
