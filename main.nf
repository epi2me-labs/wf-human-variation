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
    configure_jbrowse; } from './modules/local/common'

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
    if (!params.fast5_dir && !params.snp && !params.sv && !params.methyl && !params.cnv && !params.str) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--fast5_dir, --snp, --sv, --cnv, --str, --methyl]" + colors.reset)
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
    if(params.snp) {
        if(!params.basecaller_cfg && !params.clair3_model_path) {
            throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile> to ensure the right Clair3 model is chosen!" + colors.reset)
        }
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
    mosdepth_stats = mosdepth.out

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

    // wf-human-snp or wf-human-str
    if (params.snp || params.str) {
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
            bam_channel,
            snp_bed,
            ref_channel,
            mosdepth_stats,
            bam_stats,
            clair3_model,
        )
        output_snp(clair_vcf[0].flatten())
    }

    // wf-human-sv
    if(params.sv) {

        results = sv(
            bam_channel,
            ref_channel,
            bed,
            mosdepth_stats,
            OPTIONAL
        )
        artifacts = results[0].flatten()
        output_sv(artifacts)
    }

    // wf-human-methyl
    if (params.methyl) {
        results = methyl(bam_channel, ref_channel)
        methyl_stats = results.modbam2bed.flatten()
    }
    else {
        methyl_stats = Channel.empty()
    }

    //wf-human-cnv
    if (params.cnv) {
        results = cnv(
          bam_channel,
          bam_stats
        )
        output_cnv(results)
        }
    
    // wf-human-str
    if(params.str) {
        // use haplotagged bam from snp() as input to str()
        bam_channel_str = clair_vcf[1]

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
            methyl_stats.flatten(),
            mapula_stats.flatten(),
            jb_conf.flatten()
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
