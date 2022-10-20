#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { snp } from './workflows/wf-human-snp'
include { output_snp } from './modules/local/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { minimap2_ubam; output_sv } from './modules/local/wf-human-sv'

include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    mosdepth;
    mapula;
    getAllChromosomesBed;
    publish_artifact;
    check_for_alignment;
    configure_jbrowse; } from './modules/local/common'

include { fast5 } from './workflows/guppy'
include { methyl; output_methyl } from './workflows/methyl'

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    can_start = true
    if (!params.snp && !params.sv && !params.methyl) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --methyl]" + colors.reset)
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


    if(params.snp) {
        if(!params.model) {
            throw new Exception(colors.red + "Clair3 --model required for --snp" + colors.reset)
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

    output_bam = false // BAM/CRAM as artifact to out_dir
    is_cram = false

    if (params.fast5_dir) {
        // Basecall fast5 input
        if (params.bam) {
            throw new Exception(colors.red + "Cannot use --fast5_dir with --bam." + colors.reset)
        }
        if (!params.guppy_cfg) {
            throw new Exception(colors.red + "You must provide a guppy profile with --guppy_cfg <profile>" + colors.reset)
        }
        if (params.guppy_cfg.contains("modbase") && params.guppy_basemod_threads == 0) {
            throw new Exception(colors.red + "modbase aware guppy config requires --guppy_basemod_threads > 0" + colors.reset)
        }
        log.warn ("* Using basecalling subworkflow!")
        log.warn ("  This subworkflow is experimental and we do not provide an official guppy image at this time.")
        guppy_out = fast5(params.fast5_dir, file(params.ref)) // TODO fix ref
        bam = guppy_out.cram
        idx = guppy_out.crai
        output_bam = true // output merged CRAM to out_dir

        // construct canonical bam channel (no need for check_bam_channel here)
        bam_channel = bam.concat(idx).buffer(size: 2)
    }
    else {
        // Otherwise handle (u)BAM/CRAM
        if (!params.bam) {
            throw new Exception(colors.red + "Missing required argument --bam." + colors.reset)
        }

        check_bam = Channel.fromPath(params.bam, checkIfExists: true)
        if (params.bam.toLowerCase().endsWith("cram")) {
            is_cram = true
            check_idx = Channel.fromPath(params.bam + ".crai", checkIfExists: true)
        }
        else {
            // Just assume BAM/BAI
            check_idx = Channel.fromPath(params.bam + ".bai", checkIfExists: true)
        }
        check_bam_channel = check_bam.concat(check_idx).buffer(size: 2)
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
        check_ref = ref.concat(ref_index).buffer(size: 2) // don't wait for cram_cache to perform check_for_alignment
        checked_alignments = check_for_alignment(check_ref, check_bam_channel)

        // non-zero exit code from realignment step will be interpreted as alignment required
        // fork BAMs into realign (for (re)alignment) and noalign subchannels
        checked_alignments.branch {
            realign: it[0] == '65'
            noalign: it[0] == '0'
        }.set{alignment_fork}

        already_aligned_bams = alignment_fork.noalign.map{ it -> tuple(it[1], it[2]) }

        // Check old ref
        int to_realign = 0
        alignment_fork.realign.subscribe onNext: { to_realign++ }, onComplete: {
            if(to_realign > 0 && is_cram && !params.old_ref) {
                log.error(colors.red + "Input CRAM requires (re)alignment. You must provide a path to the original reference with '--old_ref <path>'" + colors.reset)
                throw new Exception("--old_ref not provided.")
            }
        }
        old_ref = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        if(is_cram && to_realign > 0) {
            old_ref = Channel.fromPath(params.old_ref, checkIfExists: true)
        }

        // call minimap on bams that require (re)alignment
        new_mapped_bams = minimap2_ubam(ref, old_ref, alignment_fork.realign.map{ it -> tuple(it[1], it[2]) })
        newly_aligned_bams = new_mapped_bams.alignment

        // mix realign and noalign forks back to canonical bam_channel with (reads, reads_idx) format
        bam_channel = already_aligned_bams.mix(newly_aligned_bams)
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

    // mosdepth and mapula
    mosdepth(bam_channel, bed, ref_channel)
    mosdepth_stats = mosdepth.out
    mapula(bam_channel, bed, ref_channel)
    mapula_stats = mapula.out

    // wf-human-methyl
    if (params.methyl) {
        results = methyl(bam_channel, ref_channel)
        artifacts = results[0].flatten()
        output_methyl(artifacts)
    }

    // wf-human-snp
    if (params.snp) {

        if(default_bed_set) {
            // wf-human-snp uses OPTIONAL_FILE for empty bed for legacy reasons
            snp_bed = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        }
        else {
            snp_bed = bed
        }

        model = Channel.fromPath(params.model, type: "dir", checkIfExists: true)

        clair_vcf = snp(
            bam_channel,
            snp_bed,
            ref_channel,
            mosdepth_stats,
            model
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

        if(!output_bam) {
            // filter files ending with .bam and .bam.bai, or .cram and .crai
            artifacts = artifacts.filter( { !(it.name =~ /bam($|.bai$)/) } ).filter( { !(it.name =~ /cram($|.crai$)/) } )
        }

        output_sv(artifacts)
    }

    jb_conf = configure_jbrowse(
        ref_channel,
        bam_channel,
        output_bam,
    )

    publish_artifact(
        mosdepth_stats.flatten() \
        .concat(mapula_stats.flatten()) \
        .concat(jb_conf.flatten()) \
        .concat(ref_channel.flatten())
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
