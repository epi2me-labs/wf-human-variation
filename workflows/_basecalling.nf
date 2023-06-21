include {
    merge_calls as merge_pass_calls;
    merge_calls as merge_fail_calls;
    merge_calls_to_fastq as merge_pass_calls_to_fastq;
    merge_calls_to_fastq as merge_fail_calls_to_fastq;
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
    // we can't set maxForks dynamically, but we can detect it might be wrong!
    if (task.executor != "local" && task.maxForks == 1) {
        log.warn "Non-local workflow execution detected but GPU tasks are currently configured to run in serial, perhaps you should be using '-profile discrete_gpus' to parallelise GPU tasks for better performance?"
    }
    """
    set +e
    source /opt/nvidia/entrypoint.d/*-gpu-driver-check.sh # runtime driver check msg
    set -e
    dorado basecaller \
        ${model_arg} . \
        ${remora_args} \
        ${basecaller_args} \
        --device ${params.cuda_device} | samtools view --no-PG -b -o ${chunk_idx}.ubam -
    """
}


process align_and_qsFilter {
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
        | samtools view -e '[qs] >= ${params.qscore_filter}' --output ${reads.baseName}.pass.cram --unoutput ${reads.baseName}.fail.cram -O CRAM --reference ${reference} -
    """
}


process qsFilter {
    label "wf_basecalling"
    input:
        path reads
    output:
        path("${reads.baseName}.pass.cram"), emit: pass
        path("${reads.baseName}.fail.cram"), emit: fail
    script:
    """
    samtools view -e '[qs] >= ${params.qscore_filter}' ${reads} --output ${reads.baseName}.pass.cram --unoutput ${reads.baseName}.fail.cram -O CRAM
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
        input_ref
        basecaller_model_name
        basecaller_model_path
        remora_model_name
        remora_model_path
        watch_path
        dorado_ext
        output_bam
    main:

        // determine output extentions
        def align_ext = "cram"
        def index_ext = "crai"
        if (output_bam) {
            align_ext = "bam"
            index_ext = "bai"
        }
        output_exts = Channel.of([align_ext, index_ext]).collect()

        if (input_ref) {
            if (params.fastq_only) {
                log.warn "Ignoring request to output FASTQ as you have provided a reference for alignment."
            }
            // create value channel of ref by calling first
            ref = Channel.fromPath(params.ref, checkIfExists: true).first()
        }
        else {
            ref = file("${projectDir}/data/OPTIONAL_FILE")
        }

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
        String stop_filename = "STOP.${workflow.sessionId}.${params.dorado_ext}" // use the sessionId so resume works
        existing_pod5_chunks = Channel
            .fromPath(input_path + "**.${params.dorado_ext}", checkIfExists: true)
        if (watch_path) {
            // watch input path for more pod5s
            if (params.read_limit) {
                log.warn "Watching ${input_path} for new ${dorado_ext} files, until ${params.read_limit} reads have been observed."
                log.warn "To stop this workflow early, create: ${input_path}/${stop_filename}"
            }
            else {
                log.warn "Watching ${input_path} for new ${dorado_ext} files indefinitely."
                log.warn "To stop this workflow create: ${input_path}/${stop_filename}"
            }
            watch_pod5_chunks = Channel
                .watchPath("$input_path/**.${params.dorado_ext}")
                .until{ it.name == stop_filename }
            pod5_chunks = existing_pod5_chunks
                .concat(watch_pod5_chunks)
                .buffer(size:params.basecaller_chunk_size, remainder:true)
                .map { tuple(chunk_idx++, it) }
        } else {
            pod5_chunks = existing_pod5_chunks
                .buffer(size:params.basecaller_chunk_size, remainder:true)
                .map { tuple(chunk_idx++, it) }
        }


        called_bams = dorado(
            pod5_chunks,
            tuple(basecaller_model_name, basecaller_model, basecaller_model_override),
            tuple(remora_model_name, remora_model, remora_model_override),
        )

        if (input_ref) {
            // make mmi for faster alignment
            mmi_ref = make_mmi(ref)

            // align, qscore_filter and sort
            crams = align_and_qsFilter(mmi_ref, ref, called_bams)
        }
        else {
            // skip alignment and just collate pass and fail
            crams = qsFilter(called_bams)
        }

        // merge passes and fails
        // we've aliased the merge_calls process to save writing some unpleasant looking flow
        // FASTQ output can only be used when there is no input_ref
        if (params.fastq_only && !input_ref) {
            pass = merge_pass_calls_to_fastq(crams.pass.collect(), "pass")
            fail = merge_fail_calls_to_fastq(crams.fail.collect(), "fail")
        }
        else {
            pass = merge_pass_calls(ref, crams.pass.collect(), "pass", output_exts)
            fail = merge_fail_calls(ref, crams.fail.collect(), "fail", output_exts)
        }

    emit:
        chunked_pass_crams = crams.pass
        pass = pass
        fail = fail

}
