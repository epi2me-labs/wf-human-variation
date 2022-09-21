#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { snp } from './workflows/wf-human-snp'
include { output_snp } from './modules/local/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { minimap2_ubam; output_sv } from './modules/local/wf-human-sv'

include { index_ref_gzi; index_ref_fai; cram_cache; decompress_ref } from './modules/local/common'
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

    output_bam = false

    if (params.fast5_dir) {
        // Basecall fast5 input
        if (params.bam || params.ubam ) {
            throw new Exception(colors.red + "Cannot use --fast5_dir with --bam or --ubam." + colors.reset)
        }
        log.warn ("* Using basecalling subworkflow!")
        log.warn ("  This subworkflow is not officially supported. Please do not open Github issues or raise support tickets for errors encountered attempting to use this workflow before release.")
        guppy_out = fast5(params.fast5_dir, file(params.ref)) // TODO fix ref
        bam = guppy_out.cram
        idx = guppy_out.crai
    }
    else {
        // Otherwise handle (u)BAM/CRAM
        if ((params.bam && params.ubam) || (!params.bam && !params.ubam)) {
            throw new Exception(colors.red + "Must provide one of --bam or --ubam." + colors.reset)
        }
        else if (params.bam) {
            bam = Channel.fromPath(params.bam, checkIfExists: true)
            if (params.bam.toLowerCase().endsWith("cram")) {
                idx = Channel.fromPath(params.bam + ".crai", checkIfExists: true)
            }
            else {
                // Just assume BAM/BAI
                idx = Channel.fromPath(params.bam + ".bai", checkIfExists: true)
            }
        }
        else if (params.ubam && can_start) {
            // Align uBAM with minimap2
            output_bam = true // output alignment BAM as artifact
            ubam = Channel.fromPath(params.ubam, checkIfExists: true)
            mapped = minimap2_ubam(ref, ubam)
            bam = mapped.cram
            idx = mapped.cram_index
        }
    }

    bed = null
    if(params.bed){
        bed = Channel.fromPath(params.bed, checkIfExists: true)
    }

    // Dummy optional file
    // TODO should be a channel?
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    if (can_start) {
        // Build ref cache for CRAM steps that do not take a reference
        ref_cache = cram_cache(ref)

        // Create ref index if required
        if (!ref_index_fp || !ref_index_fp.exists()) {
            log.warn ("Creating missing reference index for ${params.ref}")
            index_ref = index_ref_fai(ref)
            ref_index = index_ref.reference_index
        }
        else {
            ref_index = Channel.of(ref_index_fp)
        }

        if (params.disable_ping == false) {
            try {
                Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
            } catch(RuntimeException e1) {
            }
        }
    }

    ref_channel = ref.concat(ref_index).concat(ref_cache).buffer(size: 3)
    bam_channel = bam.concat(idx).buffer(size: 2)

    // wf-human-methyl
    if (params.methyl && can_start ) {
        results = methyl(bam_channel, ref_channel)
        artifacts = results[0].flatten()
        output_methyl(artifacts)
    }

    // wf-human-snp
    if (params.snp && can_start) {

        if(bed == null) {
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
            model
        )
        output_snp(clair_vcf[0].flatten())
    }

    // wf-human-sv
    if(params.sv && can_start) {

        results = sv(bam, idx, ref_channel, bed, OPTIONAL)
        artifacts = results[0].flatten()

        if(!output_bam) {
            // filter files ending with .bam and .bam.bai, or .cram and .crai
            artifacts = artifacts.filter( { !(it.name =~ /bam($|.bai$)/) } ).filter( { !(it.name =~ /cram($|.crai$)/) } )
        }

        output_sv(artifacts)
    }
}


if (params.disable_ping == false) {
    workflow.onComplete {
        try{
            Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }

    workflow.onError {
        try{
            Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
}
