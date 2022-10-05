
process guppy {
    label "wf_guppy"
    cpus { params.guppy_map_threads + params.guppy_basemod_threads }
    // maxForks params.gpu_count // only limit on local
    // accelerator configured by nextflow.config
    input:
        tuple val(chunk_idx), path('*')
        path ref
    output:
        tuple val(chunk_idx), path("guppy_out/pass/*.${chunk_idx}.bam"), path(ref)

    script:
    def guppy_args = params.guppy_args ?: ''
    def aligner_args = params.guppy_no_align ? '' : "-a ${ref}" // RB4F testing guppy perf on alignment
    def num_base_mod_threads = params.guppy_basemod_threads > 0 ? "--num_base_mod_threads ${params.guppy_basemod_threads}" : ''
    """
    guppy_basecaller \
        ${aligner_args} \
        -q 0 \
        --device ${params.guppy_device} \
        --num_alignment_threads ${params.guppy_map_threads} \
        ${num_base_mod_threads} \
        --config ${params.guppy_cfg} \
        --bam_out \
        ${guppy_args} \
        -i . -s guppy_out
    for bam in guppy_out/pass/*.bam ; do
        mv -v \${bam} \${bam%.bam}.${chunk_idx}.bam
    done
    """
}

process sort_calls {
    label "wf_human_sv"
    input:
        tuple val(chunk_idx), path(bams), path(ref)
    output:
        path "${chunk_idx}.cram", emit: cram
    script:
    """
    samtools cat ${bams} | samtools sort - -o ${chunk_idx}.cram -O CRAM --reference ${ref} --write-index
    """
}

process merge_calls {
    label "wf_human_sv"
    input:
        path(ref)
        path(crams)
    output:
        path "merged.cram", emit: cram
        path "merged.cram.crai", emit: crai
    script:
    """
    samtools merge merged.cram ${crams} --no-PG -O CRAM --write-index --reference ${ref}
    """
}

workflow fast5 {
    take:
        input_path
        ref
    main:

        Integer chunk_idx = 0
        def fast5_chunks = Channel
            .fromPath(input_path + "**/*.fast5")
            .buffer(size:params.fast5_chunk_size, remainder:true)
            .map { tuple(chunk_idx++, it) }

        aligned_bams = guppy(fast5_chunks, ref)
        aligned_crams = sort_calls(aligned_bams)
        out = merge_calls(ref, aligned_crams.collect())
    emit:
        cram = out.cram
        crai = out.crai
}
