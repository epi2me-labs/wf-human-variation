#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { snp } from './workflows/wf-human-snp'
include { lookup_clair3_model; output_snp } from './modules/local/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { output_sv } from './modules/local/wf-human-sv'

include { cnv } from './workflows/wf-human-cnv'
include { output_cnv } from './modules/local/wf-human-cnv'

include { output_str } from './modules/local/wf-human-str'

include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    mosdepth as mosdepth_input;
    mosdepth as mosdepth_downsampled;
    readStats;
    mapula;
    getAllChromosomesBed;
    publish_artifact;
    configure_jbrowse;
    get_coverage; 
    failedQCReport; 
    makeAlignmentReport; 
    getParams; 
    getVersions;
    getGenome; 
    eval_downsampling;
    downsampling;
    } from './modules/local/common'

include {
    bam_ingress;
} from './lib/bamingress'

include { basecalling } from './workflows/basecalling'
include { methyl; validate_modbam} from './workflows/methyl'

include { str } from './workflows/wf-human-str'

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    can_start = true
    if (!params.snp && !params.sv && !params.methyl && !params.cnv && !params.str && !params.phase_methyl) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --cnv, --str, --methyl, --phase_methyl]" + colors.reset)
        can_start = false
    }

    if (params.containsKey("ubam_threads")) {
        log.error (colors.red + "--ubam_threads is deprecated. Use `nextflow run ${workflow.manifest.name} --help` to see the parameter list." + colors.reset)
        can_start = false
    }

    if (params.containsKey("ubam")) {
        log.error (colors.red + "--ubam is deprecated as this workflow can determine whether (re)alignment is required automatically, use --bam instead." + colors.reset)
        can_start = false
    }

    // check snp has basecaller config for clair3 model lookup
    if(params.snp || params.phase_methyl) {
        if(!params.basecaller_cfg && !params.clair3_model_path) {
            throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile> to ensure the right Clair3 model is chosen!" + colors.reset)
        }
    }

    // check SV calling will be done when benchmarking SV calls
    if(params.sv_benchmark && !params.sv) {
        throw new Exception(colors.red + "Cannot benchmark SV subworkflow without running SV subworkflow! Enable the SV subworkflow with --sv." + colors.reset)
    }

    // switch workflow to BAM if calling CNV
    if (params.cnv) {
        log.warn "CNV calling subworkflow does not support CRAM. You don't need to do anything, but we're just letting you know that:"
        log.warn "- If your input file is CRAM, it will be converted to a temporary BAM inside the workflow automatically."
        log.warn "- If your input requires alignment or basecalling, the outputs will be saved to your output directory as BAM instead of CRAM."
        output_bam = true
    }
    else {
        output_bam = false
    }

    // Check ref and decompress if needed
    ref = null
    ref_index_fp = null
    if (params.ref.toLowerCase().endsWith("gz")) {
        // gzipped ref not supported by some downstream tools (pyfaidx, cram_cache)
        // easier to just decompress and pass it around rather than confusing the user
        decompress_ref(file(params.ref))
        ref = decompress_ref.out.decompressed_ref
    }
    else {
        ref = Channel.fromPath(params.ref, checkIfExists: true)
        ref_index_fp = file(params.ref + ".fai")
    }
    if (params.fast5_dir) {

        // Basecall fast5 input
        if (params.bam) {
            throw new Exception(colors.red + "Cannot use --fast5_dir with --bam." + colors.reset)
        }
        // Ensure a valid basecaller is set
        if (params.basecaller != "dorado"){
            throw new Exception(colors.red + "Basecaller ${params.basecaller} is not supported. Use --basecaller dorado." + colors.reset)
        }
        // Ensure basecaller config is set
        if (!params.basecaller_cfg && !params.basecaller_model_path) {
            throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile>" + colors.reset)
        }
        // Ensure basecaller config is not Clair3 only
        if (params.basecaller_cfg.startsWith("clair3:")) {
            throw new Exception(colors.red + "You have chosen a --basecaller_cfg that can only be used for Clair3 SNP calling.\nPlease review the list of available models with --help and pick a model that is not prefixed with 'clair3:' to enable basecalling with Dorado." + colors.reset)
        }
        if (params.basecaller_cfg == "custom" && !params.basecaller_model_path){
            throw new Exception(colors.red + "You have selected a custom basecalling model but have not provided the path of the custom model with --basecaller_model_path" + colors.reset)
        }
        if (params.remora_cfg == "custom" && !params.remora_model_path){
            throw new Exception(colors.red + "You have selected a custom modbasecalling model but have not provided the path of the custom model with --remora_model_path" + colors.reset)
        }
        // Ensure modbase threads are set if calling them
        if ((params.remora_cfg || params.remora_model_path) && params.basecaller_basemod_threads == 0) {
            throw new Exception(colors.red + "--remora_cfg modbase aware config requires setting --basecaller_basemod_threads > 0" + colors.reset)
        }

        // ring ring it's for you
        crams = basecalling(params.fast5_dir, file(params.ref), output_bam) // TODO fix file calls
        bam_channel = crams.pass
        bam_fail = crams.fail

    }
    else {
        // Set-up any non basecalling structs
        bam_fail = Channel.empty()

        // Otherwise handle (u)BAM/CRAM
        if (!params.bam) {
            throw new Exception(colors.red + "Missing required input argument, use --bam or --fast5_dir." + colors.reset)
        }

        check_bam = params.bam
    }

    // ************************************************************************
    // Bail from the workflow for a reason we should have already specified
    if (!can_start){
        throw new Exception("The workflow could not be started.")
    }
    // ************************************************************************

    // Dummy optional file
    // TODO should be a channel?
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    // Create ref index if required
    if (!ref_index_fp || !ref_index_fp.exists()) {
        index_ref = index_ref_fai(ref)
        ref_index = index_ref.reference_index
    }
    else {
        ref_index = Channel.of(ref_index_fp)
    }

    if (!params.disable_ping) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    // Determine if (re)alignment is required for input BAM
    if(!params.fast5_dir){
        bam_channel = bam_ingress(
            ref,
            ref_index,
            check_bam,
            [output_bam: output_bam],
        )
    }

    // Check if the genome build in the BAM is suitable for any workflows that have restrictions
    // NOTE getGenome will exit non-zero if the build is neither hg19 or hg38, so it shouldn't be called
    // if annotation is skipped for snp, sv and phase_methyl, to allow other genomes (including non-human)
    // to be processed

    // always getGenome for CNV and STR
    if (params.cnv || params.str) {
        genome_build = getGenome(bam_channel)
    }
    // getGenome for STP, SV and phase_methyl as long as annotation not disabled
    else if ((params.snp || params.sv || params.phase_methyl) && params.annotation) {
        genome_build = getGenome(bam_channel)
    }
    else {
        genome_build = null
    }

    // Build ref cache for CRAM steps that do not take a reference
    cram_cache(ref)
    ref_cache = cram_cache.out.ref_cache
    ref_path = cram_cache.out.ref_path
    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref.concat(ref_index).concat(ref_cache).concat(ref_path).buffer(size: 4)

    // Set BED (and create the default all chrom BED if necessary)
    bed = null
    default_bed_set = false
    if(params.bed){
        bed = Channel.fromPath(params.bed, checkIfExists: true)
    }
    else {
        default_bed_set = true
        bed = getAllChromosomesBed(ref_channel).all_chromosomes_bed
    }

    // mosdepth for depth traces -- passed into wf-snp :/
    mosdepth_input(bam_channel, bed, ref_channel)
    mosdepth_stats = mosdepth_input.out.mosdepth_tuple
    mosdepth_summary = mosdepth_input.out.summary
    if (params.depth_intervals){
        mosdepth_perbase = mosdepth_input.out.perbase
    } else {
        mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
    }
    
    // readStats for alignment and QC -- passed into wf-snp :/
    readStats(bam_channel, bed, ref_channel)
    bam_stats = readStats.out.read_stats
    bam_flag = readStats.out.flagstat

    if (params.mapula) {
        mapula(bam_channel, bed, ref_channel)
        mapula_stats = mapula.out
    }
    else {
        mapula_stats = Channel.empty()
    }
    
    // Determine if the coverage threshold is met to perform analysis.
    // If too low, it creates an empty input channel, 
    // avoiding the subsequent processes to do anything
    software_versions = getVersions()
    workflow_params = getParams()
    if (params.bam_min_coverage > 0){
        // Define if a dataset passes or not the filtering
        get_coverage(mosdepth_input.out.summary)
        // Combine with the bam and branch by passing the depth filter
        get_coverage.out.pass
            .combine(bam_channel)
            .branch{ 
                pass: it[0] == "true" 
                not_pass: it[0] == "false" 
                }
            .set{bamdepth_filter}
        // Create the pass_bam_channel  channel when they pass
        bamdepth_filter.pass
            .map{it ->
                it.size > 0 ? [it[2], it[3], it[4]] : it
            }
            .set{pass_bam_channel}
        // If it doesn't pass the minimum depth required, 
        // emit a bam channel of discarded bam files.
        bamdepth_filter.not_pass
            .subscribe {
                log.error "ERROR: File ${it[2].getName()} will not be processed by the workflow as the detected coverage of ${it[1]}x is below the minimum coverage threshold of ${params.bam_min_coverage}x required for analysis."
            }
        bamdepth_filter.not_pass
            .map{it ->
                it.size > 0 ? [it[2], it[3], it[4]] : it
            }
            .set{discarded_bams}
        discarded_bams.view()
    } else {
        // If the bam_min_depth is 0, then create alignment report for everything.
        bam_channel.set{pass_bam_channel}
        discarded_bams = Channel.empty()
    }
    // Set extensions for analyses
    extensions = pass_bam_channel.map{
        xam, xai, meta -> 
        meta.is_cram ? ['cram', 'crai'] : ['bam', 'bai']
    }

    // Check and perform downsampling if needed.
    if (params.downsample_coverage){
        // Define reduction rate
        eval_downsampling(mosdepth_input.out.summary)
        eval_downsampling.out.downsampling_ratio
            .splitCsv()
            .branch{
                subset: it[0] == 'true'
                ready: it[0] == 'false'
            }
            .set{ratio}

        // Perform downsampling
        downsampling(pass_bam_channel, ratio.subset)

        // prepare ready files
        ratio.ready
            .combine(pass_bam_channel)
            .map{ready, ratio, xam, xai, meta -> [xam, xai, meta]}
            .set{ready_bam_channel}

        // Join allowing a remainder, so that only one for each is retained.
        // we drop all null, and due to the structure the joined channel can only be:
        // - [meta, null, xam, xai], or
        // - [meta, xam, xai, null]
        // Using it - null removes the inputs from the wrong channel, retaining 
        // Before merging properly, we first check that the merged channel size is not malformed
        downsampling.out.bam
            .join(ready_bam_channel, by:2, remainder: true)
            .filter{it.size() > 4}
            .subscribe{
                throw new Exception(colors.red + "Unexpected channel size when merging." + colors.reset) 
            }
        // If this passes, then we can create the proper channel.
        downsampling.out.bam
            .join(ready_bam_channel, by:2, remainder: true)
            .map{it - null}
            .map{meta, xam, xai -> [xam, xai, meta]}
            .set{pass_bam_channel}

        // Prepare the output files for mosdepth.
        // First, we compute the depth for the downsampled files, if it
        // exists 
        mosdepth_downsampled(downsampling.out, bed, ref_channel)
        // Then, choose which output will be used in the report. 
        // If it needs to be subset, then the combined output exists, whereas 
        // the original mosdepth file is merged with the empty ready channel, leaving 
        // the correct file to output. Otherwise, the reverse happens and it emits 
        // the original mosdepth files. 
        mosdepth_summary = 
            mosdepth_downsampled.out.summary
                .combine(ratio.subset)
                .map{it[0]}
                .join(
                    mosdepth_input.out.summary
                        .combine(ratio.ready)
                        .map{it[0]}
                    , remainder: true
                    )
        mosdepth_stats = 
            mosdepth_downsampled.out.mosdepth_tuple
                .combine(ratio.subset)
                .map{[it[0], it[1], it[2]]}
                .join(
                    mosdepth_input.out.mosdepth_tuple
                        .combine(ratio.ready)
                        .map{[it[0], it[1], it[2]]}
                    , remainder: true
                    )
                .map{it - null}
        if (params.depth_intervals){
            mosdepth_perbase = 
                mosdepth_downsampled.out.perbase
                    .combine(ratio.subset)
                    .map{it[0]}
                    .join(
                        mosdepth_input.out.perbase
                            .combine(ratio.ready)
                            .map{it[0]}
                        , remainder: true
                        )
                    .map{it - null}
        } else {
            mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
        }
    }

    // Create reports for pass and fail channels
    // Create passing bam report
    report_pass = pass_bam_channel
                .combine(bam_stats)
                .combine(bam_flag)
                .combine(mosdepth_stats.map{it[0]})
                .combine(ref_channel)
                .combine(software_versions.collect())
                .combine(workflow_params)
                .flatten()
                .collect() | makeAlignmentReport
    // Create failing bam report
    report_fail = discarded_bams
                .combine(bam_stats)
                .combine(bam_flag)
                .combine(mosdepth_stats.map{it[0]})
                .combine(ref_channel)
                .combine(software_versions.collect())
                .combine(workflow_params)
                .flatten()
                .collect() | failedQCReport
    
    // wf-human-sv
    if(params.sv) {
        results = sv(
            pass_bam_channel ,
            ref_channel,
            bed,
            mosdepth_input.out.summary,
            OPTIONAL,
            genome_build,
            extensions
        )
        artifacts = results.report.flatten()
        sniffles_vcf = results.sniffles_vcf
        output_sv(artifacts)
    } else {
        sniffles_vcf = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
    }
    
    // Set up BED for wf-human-snp, wf-human-str or --phase_methyl
    if (params.snp || params.str || params.phase_methyl) {
        if(default_bed_set) {
            // wf-human-snp uses OPTIONAL_FILE for empty bed for legacy reasons
            snp_bed = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        }
        else {
            snp_bed = bed
        }

        if(params.clair3_model_path) {
            log.warn "Overriding Clair3 model with ${params.clair3_model_path}."
            clair3_model = Channel.fromPath(params.clair3_model_path, type: "dir", checkIfExists: true)
        }
        else {
            // map basecalling model to clair3 model
            lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
            // TODO basecaller_model_path
            clair3_model = lookup_clair3_model(lookup_table, params.basecaller_cfg - "clair3:").map {
                log.info "Autoselected Clair3 model: ${it[0]}" // use model name for log message
                it[1] // then just return the path to match the interface above
            }
        }

        clair_vcf = snp(
            pass_bam_channel,
            snp_bed,
            ref_channel,
            clair3_model,
            sniffles_vcf,
            genome_build,
            extensions
        )
        output_snp(clair_vcf.clair3_results.flatten())
    }

    // wf-human-methyl
    // Validate modified bam
    if (params.methyl){
        if (params.phase_methyl){
            validate_modbam(clair_vcf.hp_bams, ref_channel)
        } else {
            validate_modbam(pass_bam_channel, ref_channel)
        }

        // Warn of input without modified base tags
        validate_modbam.out.branch{
            stdbam: it[-1] == '65'
            modbam: it[-1] == '0'
            }.set{validation_results}
        // Log warn if it is not modbam
        validation_results.stdbam.subscribe{
            it -> log.warn "Input ${it[0]} does not contain modified base tags. Was a modified basecalling model selected when basecalling this data?"
        }

        // Save the other as input, keeping only the necessary elements
        validated_bam = validation_results.modbam.map{cram, crai, meta, code -> [cram, crai, meta]}
        validation_results.modbam.ifEmpty{it -> log.warn "Alignment files do not contain modified base tags. Skipping methyl aggregation and reporting."}

        results = methyl(validated_bam, ref_channel)
        methyl_stats = results.modbam2bed.flatten()
    } else {
        methyl_stats = Channel.empty()
    }

    //wf-human-cnv
    if (params.cnv) {
        results = cnv(
            pass_bam_channel,
            bam_stats,
            genome_build
        )
        output_cnv(results)
    }
    
    // wf-human-str
    if(params.str) {
        // use haplotagged bam from snp() as input to str()
        bam_channel_str = clair_vcf.str_bams

        results = str(
          bam_channel_str,
          ref_channel,
          bam_stats
        )
        output_str(results)
    }

    jb_conf = configure_jbrowse(
        ref_channel,
        bam_channel,
    )

    publish_artifact(
        // CW-1033: remove environment variable from output
        ref_channel.map{it[0..2]}.flatten().mix(
            // emit bams with the "output" meta tag
            bam_channel.filter( { it[2].output } ),
            // bam_fail can only exist if basecalling was performed
            bam_fail.flatten(),
            bam_stats.flatten(),
            mosdepth_stats.flatten(),
            mosdepth_summary.flatten(),
            mosdepth_perbase.flatten(),
            methyl_stats.flatten(),
            mapula_stats.flatten(),
            jb_conf.flatten(),
            report_pass.flatten(),
            report_fail.flatten()
        )
    )

}

if (!params.disable_ping) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
