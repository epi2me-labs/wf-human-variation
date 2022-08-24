#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { clair3 } from './workflows/wf-human-snp'
include { output_snp } from './modules/local/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { minimap2_ubam; output_sv } from './modules/local/wf-human-sv'


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    can_start = true
    if (!params.snp && !params.sv) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv]" + colors.reset)
        can_start = false
    }

    if (params.containsKey("ubam_threads")) {
        log.error (colors.red + "--ubam_threads is deprecated. Use `nextflow run ${workflow.manifest.name} --help` to see the parameter list." + colors.reset)
        can_start = false
    }

    // Check common files
    ref = Channel.fromPath(params.ref, checkIfExists: true)

    output_bam = false
    if ((params.bam && params.ubam) || (!params.bam && !params.ubam)) {
        throw new Exception("Must provide one of --bam or --ubam.")
        can_start = false
    }
    else if (params.bam) {
        bam = Channel.fromPath(params.bam, checkIfExists: true)
        bai = Channel.fromPath(params.bam + ".bai", checkIfExists: true)
    }
    else if (params.ubam && can_start) {
        output_bam = true // output alignment BAM as artifact
        ubam = Channel.fromPath(params.ubam, checkIfExists: true)
        mapped = minimap2_ubam(ref, ubam)
        bam = mapped.bam
        bai = mapped.bam_index
    }


    bed = null
    if(params.bed){
        bed = Channel.fromPath(params.bed, checkIfExists: true)
    }

    // Dummy optional file
    // TODO should be a channel?
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    if (can_start) {
        if (params.disable_ping == false) {
            try {
                Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
            } catch(RuntimeException e1) {
            }
        }
    }

    if (params.snp && can_start) {
        fai = Channel.fromPath(params.ref + ".fai", checkIfExists: true)
        snp_ref = ref.concat(fai).buffer(size: 2)
        snp_bam = bam.concat(bai).buffer(size: 2)

        if(bed == null) {
            snp_bed = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        }
        else {
            snp_bed = bed
        }
        model = Channel.fromPath(params.model, type: "dir", checkIfExists: true)

        clair_vcf = clair3(
            snp_bam,
            snp_bed,
            snp_ref,
            model
        )
        output_snp(clair_vcf[0].flatten())
    }
    if(params.sv && can_start) {

        results = sv(bam, bai, ref, bed, OPTIONAL)
        artifacts = results[0].flatten()

        if(!output_bam) {
            // filter files ending with .bam and .bam.bai
            artifacts = artifacts.filter( { !(it.name =~ /bam($|.bai$)/) } )
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
