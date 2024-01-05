process modkit {
    label "wf_human_mod"
    cpus params.modkit_threads
    input:
        tuple path(alignment), path(alignment_index), val(alignment_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    output:
        path "${params.sample_name}*.bedmethyl.gz", emit: modkit_outputs

    script:
    """
    modkit pileup \\
        ${alignment} \\
        ${params.sample_name}.wf_mods.bedmethyl \\
        --ref ${ref} \\
        --threads ${task.cpus} ${options}
    
    # Compress all
    bgzip ${params.sample_name}.wf_mods.bedmethyl
    """
}

process modkit_phase {
    label "wf_human_mod"
    cpus params.modkit_threads
    input:
        tuple path(alignment), path(alignment_index), val(alignment_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    output:
        path "${params.sample_name}/${params.sample_name}*.bedmethyl.gz", emit: modkit_outputs

    script:
    // CW-2370: modkit saves in a directory when using --partition-tag rather than a single file
    // Post processing fixes this by compressing the different files and moving them to the correct directory.
    """
    modkit pileup \\
        ${alignment} \\
        ${params.sample_name} \\
        --ref ${ref} \\
        --partition-tag HP \\
        --prefix ${params.sample_name} \\
        --threads ${task.cpus} ${options}
    
    # Compress all
    for i in `ls ${params.sample_name}/`; do
        root_name=\$( basename \$i '.bed' )
        mv ${params.sample_name}/\${root_name}.bed ${params.sample_name}/\${root_name}.wf_mods.bedmethyl
        bgzip ${params.sample_name}/\${root_name}.wf_mods.bedmethyl
    done
    """
}

// Check that the bam has modifications
process validate_modbam {
    input:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta) 
        tuple path(reference), 
            path(reference_index), 
            path(reference_cache),
            env(REF_PATH)
    output:
        tuple path(alignment), 
            path(alignment_index), 
            val(meta),
            env(valid)

    script:
    """
    valid=0
    workflow-glue check_valid_modbam ${alignment} || valid=\$?

    # Allow EX_OK and EX_DATAERR, otherwise explode
    if [ \$valid -ne 0 ] && [ \$valid -ne 65 ]; then
        exit 1
    fi
    """
}


workflow mod {
    take:
        alignment
        reference
    main:
        def modkit_options = params.force_strand ? '' : '--combine-strands --cpg'
        // Custom options overwrite every custom setting
        if (params.modkit_args){
            modkit_options = "${params.modkit_args}"
        }

        // CW-2370: modkit doesn't require to treat each haplotype separately, as
        // you simply provide --partition-tag HP and it will automatically generate
        // three distinct output files, one for each haplotype and one for the untagged regions.
        if (params.phased){
            out = modkit_phase(alignment, reference.collect(), modkit_options)
        } else {
            out = modkit(alignment, reference.collect(), modkit_options)
        }
    emit:
        modkit = out
}
