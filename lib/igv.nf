/*
 This workflow adds the IGV configuration for
 the wf-*-variation family.
*/
include { configure_igv } from './common'
workflow igv {
    take:
        outputs
    main:
        // Create IGV configuration
        // To avoid saving the reference as an output (we just recently removed it),
        // We instead pass it as a string of params.ref. Given the following assumptions:
        // 1. the igv.json is only relevant in the desktop app, and
        // 2. the desktop app only passes absolute paths when referring to files outside the outdir
        // The IGV configuration should work even when passing the input variable as file name.
        boolean keep_track_order = false
        igv_files = outputs
            | flatten
            // If a file path starts with the workdir, then we assume it will be output by the workflow in the outdir and refer to it by its file name
            // Otherwise, we assume the path is an existing input and refer to it with the file path.
            | map {
                it ->
                def outfile = it.toString().startsWith("${workflow.workDir}") ? "${it.name}" : "${it.toString()}"
                outfile
            }
            | collectFile(name: "file-names.txt", newLine: true, sort: false)

        // Create IGV configuration file.
        igv_out = configure_igv(igv_files, Channel.of(null), Channel.of(null), Channel.of(null), keep_track_order)
    emit:
        igv_out = igv_out
}