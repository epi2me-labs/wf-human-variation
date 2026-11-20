// Pre-compute sample probabilities
process sample_probs {
    label "wf_human_mod"
    // Using 4 threads on a 90X takes ~30sec to complete
    cpus 4
    memory { 15.GB * task.attempt - 1.GB }
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
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
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta), path(xam), path(xai)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    output:
        tuple val(meta), val('*'), path("${meta.alias}*bedmethyl.gz"), emit: modkit

    script:
    """
    modkit pileup \\
        ${xam} \\
        ${meta.alias}.wf_mods.${meta.sq}.bedmethyl.gz \\
        --modified-bases 6mA \\
        --ref ${ref} \\
        --region ${meta.sq} \\
        --log-filepath modkit.log \\
        ${meta.probs} \\
        --bgzf \\
        --threads ${task.cpus} ${options}
    """
}

process modkit_phase {
    label "wf_human_mod"
    cpus params.modkit_threads
    // Phasing is a bit more greedy for memory. Use 2.GB/core + buffer.
    memory {(2.GB * params.modkit_threads * task.attempt) + 3.GB}
    maxRetries 1
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple val(meta), path(xam), path(xai)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val options
    // some of the outputs can be optional based on the tagging (they can all be from one hap, either haps, or none)
    output:
        tuple val(meta), val('*'), path("${meta.alias}/${meta.alias}*combined.bedmethyl.gz"), emit: modkit_Hstar, optional: true
        tuple val(meta), val('1'), path("${meta.alias}/${meta.alias}*1.bedmethyl.gz"), emit: modkit_H1, optional: true
        tuple val(meta), val('2'), path("${meta.alias}/${meta.alias}*2.bedmethyl.gz"), emit: modkit_H2, optional: true

    script:
    // CW-2370: modkit saves in a directory when using --partition-tag rather than a single file
    // Post processing fixes this by compressing the different files and moving them to the correct directory.
    """
    modkit pileup \\
        ${xam} \\
        ${meta.alias} \\
        --modified-bases 6mA \\
        --ref ${ref} \\
        --phased \\
        --prefix ${meta.alias}.wf_mods.${meta.sq} \\
        --log-filepath modkit.log \\
        --region ${meta.sq} \\
        ${meta.probs} \\
        --bgzf \\
        --threads ${task.cpus} ${options}
    
    # Compress all
    for i in `ls ${meta.alias}/`; do
        root_name=\$( basename \$i '.bed.gz' )
        # modkit saves the file as meta.alias.wf_mods.sq_haplotype.bed.gz
        # create a new name with the patter meta.alias.wf_mods.haplotype.bedmethyl
        new_name=\$( echo \${root_name} | sed 's/wf_mods\\.${meta.sq}_/wf_mods\\.${meta.sq}\\./' )
        mv ${meta.alias}/\${root_name}.bed.gz ${meta.alias}/\${new_name}.bedmethyl.gz
    done
    """
}

process concat_bedmethyl {
    cpus 4
    memory 8.GB

    input:
        tuple val(meta), val(group), path("bedmethyls/*")
    output:
        tuple val(meta), val(group), path("${meta.alias}.wf_mods.*bedmethyl.gz")

    script:
    // Concatenate the bedMethyl, sort them and compress them
    def label = group != '*' ? "${group}." : ""
    """
    zcat -f bedmethyls/* | \
        sort -k 1,1 -k2,2n --parallel ${task.cpus} | \
        bgzip -c -@ ${task.cpus} > ${meta.alias}.wf_mods.${label}bedmethyl.gz
    """
}

process modkit_tobigwig {
    label "wf_human_mod"
    cpus 4
    memory 2.GB
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(meta), val(group), path(bedmethyl), val(mod)
    output:
        tuple val(meta), val(group), val(mod), path("${meta.alias}.wf_mods.*.bw")
    publishDir \
        path: "${params.out_dir}",
        mode: 'copy'

    script:
    def mod_code
    def mod_label
    if(mod.contains(':')) {
        (mod_code, mod_label) = mod.split(':', 2)
    }
    else {
        mod_code = mod_label = mod
    }
    def label = group != '*' ? "${group}-${mod_label}" : "${mod_label}"
    // switch on inverted counting for negative strand if using force_strand
    def strand_values_arg = params.force_strand ? "--negative-strand-values" : ""
    """
    zcat ${bedmethyl} | modkit bm tobigwig --sizes ${ref_idx} -t ${task.cpus} --mod-codes ${mod_code} ${strand_values_arg} - ${meta.alias}.wf_mods.${label}.bw
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
        run_haplotagging
    main:
        def modkit_options = ''
        // Custom options overwrite every custom setting.
        if (params.modkit_args){
            modkit_options = "${params.modkit_args}"
        }

        // Create input channel
        // Process only contigs with reads mapped to them.
        // Follow behaviour of SNP and SV by filtering contigs to chromosomes_codes,
        // if --`include_all_ctgs false`.
        // Add chromosome and filtering probs to the modkit channel.
        // If haplotagging is requested, then the channel already has the chromosome ID
        // in it, so ignore the target_chrom.
        if (run_haplotagging){
            modkit_bam = modkit_bam
            | combine(probs)
            | map{
                xam, xai, meta, probs ->
                [meta + [probs: probs], xam, xai]
            }
        } else {
            target_chrom = bam_flagstats
                | splitCsv(sep: "\t", header: true)
                | filter{
                    it.ref != "*" && it.total as int > 0
                }
                | filter{
                    params.include_all_ctgs ? true : chromosome_codes.contains(it.ref)
                }
                | map{it.ref}
            modkit_bam = target_chrom
            | combine(modkit_bam)
            | combine(probs)
            | map{
                sq, xam, xai, meta, probs ->
                [meta + [sq:sq, probs:probs], xam, xai]
            }
        }

        // CW-2370: modkit doesn't require to treat each haplotype separately, as
        // you simply provide --partition-tag HP and it will automatically generate
        // three distinct output files, one for each haplotype and one for the untagged regions.
        if (params.phased){
            // Process the chunked haplotagged BAM file.
            modkit_out = modkit_phase(modkit_bam, reference.collect(), modkit_options)
            // Concatenate the haplotypes.
            bedmethyl = modkit_out.modkit_Hstar
                | mix(modkit_out.modkit_H1, modkit_out.modkit_H2)
                | map{ meta, group, bedmethyl -> [["alias": meta.alias], group, bedmethyl]}
                | groupTuple(by: [0,1])
                | concat_bedmethyl
        } else {
            // Run modkit.
            bedmethyl = modkit(modkit_bam, reference.collect(), modkit_options)
                | map{ meta, group, bedmethyl -> [["alias": meta.alias], group, bedmethyl]}
                | groupTuple(by: [0,1])
                | concat_bedmethyl
        }
        // make bigwig from the (un)phased bedmethyls
        // Generate a 6mA bigWig track for adenine modifications.
        // mod codes use basemod codes from SAMTags, we label the files with the more
        //  human friendly abbreviation by sending <mod_code>:<mod_label> to tobigwig
        bigwig = modkit_tobigwig(
            reference.collect(),
            bedmethyl | combine(["a:6mA"]),
        )

        // IGV tracks
        // send H0, unless phased then send H1, H2
        // TODO(CW-6944) send bigwig files to IGV instead
        igv = bedmethyl
          | filter { _meta, group, _bedmethyl ->
              (!params.phased && group == "*") || (params.phased && group != "*")
          }
          | map { _meta, _group, bedmethyl -> bedmethyl }
          | collect
    emit:
        bedmethyl = bedmethyl | map { _meta, _group, bedmethyl -> bedmethyl } | collect
        bigwig = bigwig | map { _meta, _group, _mod, bigwig -> bigwig } | collect
        igv = igv
}
