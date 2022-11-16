
process dorado {
    label "wf_dorado"
    label "wf_basecalling"
    label "gpu"
    accelerator 1 // further configuration should be overloaded using withLabel:gpu
    cpus 8
    input:
        tuple val(chunk_idx), path('*')
        path ref
    output:
        path("${chunk_idx}.ubam")
    script:
    def remora_args = (params.basecaller_basemod_threads > 0 && params.remora_cfg) ? "--remora-models \${DRD_MODELS_PATH}/${params.remora_cfg} --remora-threads ${params.basecaller_basemod_threads} --remora-batchsize 1024" : ''
    """
    dorado basecaller \
        "\${DRD_MODELS_PATH}/${params.basecaller_cfg}" . \
        ${remora_args} \
        --device ${params.cuda_device} | samtools view -b -o ${chunk_idx}.ubam -
    """
}


process dorado_align {
    label "wf_basecalling"
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        path mmi_reference
        path reference
        path reads
    output:
        // NOTE merge does not need an index if merging with region/BED (https://github.com/samtools/samtools/blob/969d44990df7fa9c7bda3a7140a2c1d1bd8c62a0/bam_sort.c#L1256-L1272)
        // so we can save a few cycles and just output the CRAM
        path("${reads.baseName}.cram")
    script:
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${reads} \
        | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${mmi_reference} - \
        | samtools sort -@ ${params.ubam_sort_threads} -o ${reads.baseName}.cram -O CRAM --reference ${reference} -
    """
}


process merge_calls {
    label "wf_basecalling"
    cpus 4
    input:
        path(ref)
        path(crams)
    output:
        path "${params.sample_name}.cram", emit: cram
        path "${params.sample_name}.cram.crai", emit: crai
    script:
    """
    samtools merge ${params.sample_name}.cram ${crams} --no-PG -O CRAM --write-index --reference ${ref} --threads ${task.cpus}
    """
}


process make_mmi {
    label "wf_basecalling"
    cpus 4
    input:
        path(ref)
    output:
        path("ref.mmi")
    script:
    """
    minimap2 -t ${task.cpus} -x map-ont -d ref.mmi ${ref}
    """
}


workflow wf_dorado {
    take:
        input_path
        ref
    main:
        Integer chunk_idx = 0
        def pod5_chunks = Channel
            .fromPath(input_path + "**.${params.dorado_ext}")
            .buffer(size:params.basecaller_chunk_size, remainder:true)
            .map { tuple(chunk_idx++, it) }
        called_bams = dorado(pod5_chunks, ref)

        // make mmi for faster alignment
        mmi_ref = make_mmi(ref)

        // align and sort
        aligned_crams = dorado_align(mmi_ref, ref, called_bams)
        out = merge_calls(ref, aligned_crams.collect())
    emit:
        cram = out.cram
        crai = out.crai
}
