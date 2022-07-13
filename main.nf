#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

nextflow.enable.dsl = 2

include { start_ping; end_ping } from './lib/ping'

include { clair3 } from './workflows/wf-human-snp'
include { output_snp } from './modules/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { minimap2_ubam; output_sv } from './modules/wf-human-sv'

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    start_ping()

    // Check common files
    ref = Channel.fromPath(params.ref, checkIfExists: true)

    output_bam = false
    if ((params.bam && params.ubam) || (!params.bam && !params.ubam)) {
        println("Must provide one of --bam or --ubam.")
        exit 1
    }
    else if (params.bam) {
        bam = Channel.fromPath(params.bam, checkIfExists: true)
        bai = Channel.fromPath(params.bam + ".bai", checkIfExists: true)
    }
    else if (params.ubam) {
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

    if (params.snp) {
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
    if(params.sv) {

        results = sv(bam, bai, ref, bed, OPTIONAL)
        artifacts = results[0].flatten()

        if(!output_bam) {
            // filter files ending with .bam and .bam.bai
            artifacts = artifacts.filter( { !(it.name =~ /bam($|.bai$)/) } )
        }

        output_sv(artifacts)
    }

    //TODO end ping telemetry
    end_ping(OPTIONAL)
}
