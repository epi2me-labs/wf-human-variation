include {
    merge_calls as merge_pass_calls;
    merge_calls as merge_fail_calls;
} from './merge.nf'  // sigh. module aliasing seems easier than flow.


process dorado {
    label "wf_dorado"
    label "wf_basecalling"
    label "gpu"
    accelerator 1 // further configuration should be overloaded using withLabel:gpu
    cpus 8
    input:
        tuple val(chunk_idx), path('*')
        tuple val(basecaller_cfg), path("dorado_model"), val(basecaller_model_override)
        tuple val(remora_cfg), path("remora_model"), val(remora_model_override)
    output:
        path("${chunk_idx}.ubam")
    script:
    def remora_model = remora_model_override ? "remora_model" : "\${DRD_MODELS_PATH}/${remora_cfg}"
    def remora_args = (params.basecaller_basemod_threads > 0 && (params.remora_cfg || remora_model_override)) ? "--modified-bases-models ${remora_model}" : ''
    def model_arg = basecaller_model_override ? "dorado_model" : "\${DRD_MODELS_PATH}/${basecaller_cfg}"
    def basecaller_args = params.basecaller_args ?: ''
    """
    echo '***'
    echo 'Available models:'
    list-models | sed 's,^,- ,' | sed "s,\${DRD_MODELS_PATH}/,,"
    echo '***'
    echo 'You selected:'
    echo "Basecalling model: ${basecaller_cfg}"
    echo "Remora model     : ${remora_cfg}"
    echo '***'
    echo 'A file open error below indicates that you have entered an unknown model name.'
    echo 'It is possible the model you selected worked previously but has been updated to a new version.'
    echo 'Resubmit this workflow with an appropriate model from the model list above.'
    echo '***'

    dorado basecaller \
        ${model_arg} . \
        ${remora_args} \
        ${basecaller_args} \
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
        path("${reads.baseName}.pass.cram"), emit: pass
        path("${reads.baseName}.fail.cram"), emit: fail
    script:
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${reads} \
        | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${mmi_reference} - \
        | samtools sort -@ ${params.ubam_sort_threads} \
        | tee >(samtools view -e '[qs] < ${params.qscore_filter}}' -o ${reads.baseName}.fail.cram -O CRAM --reference ${reference} - ) \
        | samtools view -e '[qs] >= ${params.qscore_filter}' -o ${reads.baseName}.pass.cram -O CRAM --reference ${reference} -
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
        basecaller_model_name
        basecaller_model_path
        remora_model_name
        remora_model_path
    main:
        // Munge models
        // I didn't want to use the same trick from wf-humvar as I thought the models here are much larger
        // ...they aren't, but nevermind this is less hilarious than the humvar way
        basecaller_model = file("${projectDir}/data/OPTIONAL_FILE")
        def basecaller_model_override = false
        if (params.basecaller_model_path) {
            basecaller_model = file(params.basecaller_model_path, type: "dir", checkIfExists: true)
            basecaller_model_override = true
            log.warn "Overriding basecaller model with ${params.basecaller_model_path}"
        }
        remora_model = file("${projectDir}/data/OPTIONAL_FILE")
        def remora_model_override = false
        if (params.remora_model_path) {
            remora_model = file(params.remora_model_path, type: "dir", checkIfExists: true)
            remora_model_override = true
            log.warn "Overriding remora model with ${params.remora_model_path}"
        }

        Integer chunk_idx = 0
        pod5_chunks = Channel
            .fromPath(input_path + "**.${params.dorado_ext}", checkIfExists: true)
            .buffer(size:params.basecaller_chunk_size, remainder:true)
            .map { tuple(chunk_idx++, it) }
        called_bams = dorado(
            pod5_chunks,
            tuple(basecaller_model_name, basecaller_model, basecaller_model_override),
            tuple(remora_model_name, remora_model, remora_model_override),
        )

        // make mmi for faster alignment
        mmi_ref = make_mmi(ref)

        // align, qscore_filter and sort
        aligned_crams = dorado_align(mmi_ref, ref, called_bams)

        // merge passes and fails
        // we've aliased the merge_calls process to save writing some
        // unpleasant looking flow
        pass = merge_pass_calls(ref, aligned_crams.pass.collect(), "pass")
        fail = merge_fail_calls(ref, aligned_crams.fail.collect(), "fail")
    emit:
        pass = pass
        fail = fail
}
