// Pre-compute sample probabilities
process sample_probs {
    label "wf_human_mod"
    // Using 4 threads on a 90X takes ~30sec to complete
    cpus 4
    memory 8.GB
    input:
        tuple path(xam), path(xam_index), val(meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        env(probs), emit: probs

    script:
    // Set `--interval-size` to 5Mb to speed up sampling, and `--only-mapped -p 0.1` to be consistent with `pileup`
    """
    probs=\$( modkit sample-probs ${xam} -p 0.1 --interval-size 5000000 --only-mapped --threads ${task.cpus} 2> /dev/null | awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' )
    """
}

process modkit {
    label "wf_human_mod"
    cpus params.modkit_threads
    memory {(1.GB * params.modkit_threads * task.attempt) + 3.GB}
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(chr), path(alignment), path(alignment_index), val(alignment_meta), val(probs)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    output:
        tuple val('*'), path("${params.sample_name}*bedmethyl.gz"), emit: modkit

    script:
    """
    modkit pileup \\
        ${alignment} \\
        ${params.sample_name}.wf_mods.${chr}.bedmethyl \\
        --ref ${ref} \\
        --region ${chr} \\
        --interval-size 1000000 \\
        --log-filepath modkit.log \\
        ${probs} \\
        --threads ${task.cpus} ${options}
    
    # Compress all
    bgzip ${params.sample_name}.wf_mods.${chr}.bedmethyl
    """
}

process modkit_phase {
    label "wf_human_mod"
    cpus params.modkit_threads
    // Phasing is a bit more greedy for memory. Use 2.GB/core + buffer.
    memory {(2.GB * params.modkit_threads * task.attempt) + 3.GB}
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(chr), path(xam), path(xam_idx), val(probs)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    // some of the outputs can be optional based on the tagging (they can all be from one hap, either haps, or none)
    output:
        tuple val('ungrouped'), path("${params.sample_name}/${params.sample_name}*ungrouped.bedmethyl.gz"), emit: modkit_H0, optional: true
        tuple val('1'), path("${params.sample_name}/${params.sample_name}*1.bedmethyl.gz"), emit: modkit_H1, optional: true
        tuple val('2'), path("${params.sample_name}/${params.sample_name}*2.bedmethyl.gz"), emit: modkit_H2, optional: true

    script:
    // CW-2370: modkit saves in a directory when using --partition-tag rather than a single file
    // Post processing fixes this by compressing the different files and moving them to the correct directory.
    """
    modkit pileup \\
        ${xam} \\
        ${params.sample_name} \\
        --ref ${ref} \\
        --partition-tag HP \\
        --interval-size 1000000 \\
        --prefix ${params.sample_name}.wf_mods.${chr} \\
        --log-filepath modkit.log \\
        --region ${chr} \\
        ${probs} \\
        --threads ${task.cpus} ${options}
    
    # Compress all
    for i in `ls ${params.sample_name}/`; do
        root_name=\$( basename \$i '.bed' )
        # modkit saves the file as params.sample_name.wf_mods_haplotype.bed
        # create a new name with the patter params.sample_name.wf_mods.haplotype.bedmethyl
        new_name=\$( echo \${root_name} | sed 's/wf_mods_/wf_mods\\./' )
        mv ${params.sample_name}/\${root_name}.bed ${params.sample_name}/\${new_name}.bedmethyl
        bgzip ${params.sample_name}/\${new_name}.bedmethyl
    done
    """
}

process concat_bedmethyl {
    cpus 4
    memory 8.GB

    input:
        tuple val(group), path("bedmethyls/*")
    output:
        path "${params.sample_name}.wf_mods.*bedmethyl.gz"

    script:
    // Concatenate the bedMethyl, sort them and compress them
    def label = group != '*' ? "${group}." : ""
    """
    zcat -f bedmethyls/* | \
        sort -k 1,1 -k2,2n --parallel ${task.cpus} | \
        bgzip -c -@ ${task.cpus} > ${params.sample_name}.wf_mods.${label}bedmethyl.gz
    """
}

// Check that the bam has modifications
process validate_modbam {
    label "wf_common"
    cpus 1
    memory 4.GB
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
        modkit_bam
        bam_flagstats
        chromosome_codes
        probs
        reference
    main:
        def modkit_options = params.force_strand ? '' : '--combine-strands --cpg'
        // Custom options overwrite every custom setting.
        if (params.modkit_args){
            modkit_options = "${params.modkit_args}"
        }

        // CW-2370: modkit doesn't require to treat each haplotype separately, as
        // you simply provide --partition-tag HP and it will automatically generate
        // three distinct output files, one for each haplotype and one for the untagged regions.
        if (params.phased){
            // Process the chunked haplotagged BAM file.
            modkit_out = modkit_phase(modkit_bam.combine(probs), reference.collect(), modkit_options)
            // Concatenate the haplotypes.
            out = modkit_out.modkit_H0 
                | mix(modkit_out.modkit_H1, modkit_out.modkit_H2)
                | groupTuple(by: 0)
                | concat_bedmethyl
        } else {
            // Process only contigs with reads mapped to them.
            // Follow behaviour of SNP and SV by filtering contigs to chromosomes_codes,
            // if --`include_all_ctgs false`.
            target_chrom = bam_flagstats
                | splitCsv(sep: "\t", header: true)
                | filter{
                    it.ref != "*" && it.total as int > 0
                }
                | filter{
                    params.include_all_ctgs ? true : chromosome_codes.contains(it.ref)
                }
                | map{it.ref}
            // Add chromosome and filtering probs to the modkit channel.
            modkit_bam = target_chrom.combine(modkit_bam).combine(probs)
            // Run modkit.
            out = modkit(modkit_bam, reference.collect(), modkit_options)
                | groupTuple(by: 0)
                | concat_bedmethyl
        }
    emit:
        modkit = out
}
