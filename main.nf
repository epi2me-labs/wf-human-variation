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
    mosdepth;
    readStats;
    mapula;
    getAllChromosomesBed;
    publish_artifact;
    configure_jbrowse;
    get_coverage; 
    failedQCReport; 
    getParams; 
    getVersions;
    getGenome; } from './modules/local/common'

include {
    bam_ingress;
} from './lib/bamingress'

include { basecalling } from './workflows/basecalling'
include { methyl; } from './workflows/methyl'

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
        // Ensure modbase threads are set if calling them
        if ((params.remora_cfg || params.remora_model_path) && params.basecaller_basemod_threads == 0) {
            throw new Exception(colors.red + "--remora_cfg modbase aware config requires setting --basecaller_basemod_threads > 0" + colors.reset)
        }

        // ring ring it's for you
        crams = basecalling(params.fast5_dir, file(params.ref)) // TODO fix file calls
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
        )
    }

    // Check if the genome build in the BAM is suitable for any workflows that have restrictions
    // NOTE getGenome will exit non-zero if the build is neither hg19 or hg38, so don't call it unless needed
    if (params.str || params.cnv){
        // CNV requires hg19 or hg38
        // STR requires hg38
        genome_build = getGenome(bam_channel)
    }
    else {
        genome_build = null // only CNV consumes genome_build currently
    }

    // Build ref cache for CRAM steps that do not take a reference
    ref_cache = cram_cache(ref)

    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref.concat(ref_index).concat(ref_cache).buffer(size: 3)

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
    mosdepth(bam_channel, bed, ref_channel)
    mosdepth_stats = mosdepth.out.mosdepth_tuple
    mosdepth_perbase = mosdepth.out.perbase

    // readStats for alignment and QC -- passed into wf-snp :/
    readStats(bam_channel, bed, ref_channel)
    bam_stats = readStats.out

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
    if (params.bam_min_coverage > 0){
        // Define if a dataset passes or not the filtering
        get_coverage(mosdepth.out.summary)
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
            .set{pass_bam_channel }
        // If it doesn't pass the minimum depth required, 
        // emit a bam channel of discarded bam files.
        bamdepth_filter.not_pass
            .map{it ->
                it.size > 0 ? [it[1], it[2], it[3], it[4]] : it
            }
            .set{discarded_bams}
        discarded_bams
            .subscribe {
                log.error "ERROR: File ${it[1].getName()} will not be processed by the workflow as the detected coverage of ${it[0]}x is below the minimum coverage threshold of ${params.bam_min_coverage}x required for analysis."
            }
        // Create a report if the discarded bam channel is not empty.
        software_versions = getVersions()
        workflow_params = getParams()
        report = failedQCReport(
            discarded_bams, bam_stats, mosdepth_stats,
            software_versions.collect(), workflow_params)
    } else {
        // If the bam_min_depth is 0, then run as usual.
        bam_channel.set{pass_bam_channel }
        report = Channel.empty()
    }

    // wf-human-sv
    if(params.sv) {
        results = sv(
            pass_bam_channel ,
            ref_channel,
            bed,
            mosdepth_stats,
            OPTIONAL
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
            clair3_model = lookup_clair3_model(lookup_table, params.basecaller_cfg)
        }

        clair_vcf = snp(
            pass_bam_channel,
            snp_bed,
            ref_channel,
            mosdepth_stats,
            bam_stats,
            clair3_model,
            sniffles_vcf
        )
        output_snp(clair_vcf.clair3_results.flatten())
    }

    // wf-human-methyl
    if (params.methyl || params.phase_methyl) {
        if (params.phase_methyl){
            results = methyl(clair_vcf.hp_bams, ref_channel)
        } else {
            results = methyl(pass_bam_channel, ref_channel)
        }
        methyl_stats = results.modbam2bed.flatten()
    }
    else {
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
        ref_channel.flatten().mix(
            // emit bams with the "output" meta tag
            bam_channel.filter( { it[2].output } ),
            // bam_fail can only exist if basecalling was performed
            bam_fail.flatten(),
            bam_stats.flatten(),
            mosdepth_stats.flatten(),
            mosdepth_perbase.flatten(),
            methyl_stats.flatten(),
            mapula_stats.flatten(),
            jb_conf.flatten(),
            report.flatten()
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
