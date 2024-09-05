// This subworkflow aims at pre-processing an input reference
// genome, preparing the appropriate indexes and cache.
// By default, the workflow takes an input reference,
// decompress it if it is gzipped, and import (if available)
// or index with `samtools faidx` the reference genome.
// if the user requests it, the workflow can also:
// 1. generate the cram cache, with annex environmental
//    `REF_PATH` variable
// 2. generate the minimap2 `.mmi` index, for faster alignment
Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

// Argument parser
Map parse_reference(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args:[
            "input_ref",
        ],
        kwargs:[
            "output_cache": false,
            "output_mmi": false,
        ],
        name: "reference_ingress")
    return parser.parse_args(arguments)
}

// Process to generate the CRAM cache and
// create the REF_PATH variable
process cram_cache {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path reference
    output:
        tuple path("ref_cache/"), env(REF_PATH), emit: ref_cache
    shell:
    '''
    # Invoke from binary installed to container PATH
    seq_cache_populate.pl -root ref_cache/ !{reference}
    REF_PATH="ref_cache/%2s/%2s/%s"
    '''
}

// Process to create the faidx index
process faidx {
    // If the input file is gzipped, we need to emit the indexes for the input gzip file
    // only. Therefore, this become reduntant to be emitted as it won't be used by the
    // IGV configuration, but only by internal processes together with the decompressed
    // FASTA file. To avoid unnecessary emissions, we enable only if the input file is
    // decompressed.
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", enabled: !params.ref.toLowerCase().endsWith("gz")
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path(ref)
    output:
        path("${ref}.fai")
    script:
    """
    samtools faidx ${ref}
    """
}

// Process to create the faidx indexes for a gzipped reference
process gz_faidx {
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    label "wf_common"
    cpus 1
    memory 4.GB
    // If a user provides a non-bgzipped file, the process won't
    // generate the indexes. We should tolerate that, still avoid emitting
    // the reference and simply have a broken IGV file.
    // The gzi is not required to operate the workflow, so we actually tolerate any failure.
    errorStrategy 'ignore'
    input:
        path(ref)
    output:
        tuple path(ref), path("${ref}.fai"), path("${ref}.gzi")
    script:
    """
    samtools faidx ${ref}
    """
}

// Decompress the reference genome, if it is compressed
// NOTE -f required to compress symlink
process decompress_ref {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path "ref.fa.gz"
    output:
        path "ref.fa", emit: decompressed_ref
    """
    gzip -df ref.fa.gz
    """
}

// Prepare minimap2 .mmi index
process make_mmi {
    // Minimap2 is not available in wf_common
    cpus 4
    memory 16.GB
    input:
        path("ref.fa")
    output:
        path("ref.mmi")
    script:
    """
    minimap2 -t ${task.cpus} -x map-ont -d ref.mmi ref.fa
    """
}


// Workflow to prepare the reference genome and its indexes.
workflow prepare_reference {
    take:
        arguments
    main:
        Map margs = parse_reference(arguments)

        // Base ref channel
        if (file(margs.input_ref).exists()){
            ref = Channel.fromPath(margs.input_ref)
        } else {
            throw new Exception(colors.red + "File ${margs.input_ref} not found." + colors.reset)
        }

        // Check ref and decompress if needed
        // gzipped ref not supported by some downstream tools (e.g. cram_cache)
        // easier to just decompress and pass it around rather than confusing the user
        is_compressed = margs.input_ref.toLowerCase().endsWith("gz")
        if (is_compressed) {
            // Define indexes names.
            String input_fai_index = "${margs.input_ref}.fai"
            String input_gzi_index = "${margs.input_ref}.gzi"

            // Decompress reference genome.
            ref = decompress_ref(ref)
            // Check whether the input gzref is indexed. If so, pass these as indexes.
            // Otherwise, generate the gzip + fai indexes for the compressed reference.
            if (file(input_fai_index).exists() && file(input_gzi_index).exists()){
                gzindexes = Channel.fromPath(margs.input_ref)
                | mix(
                    Channel.fromPath(input_fai_index),
                    Channel.fromPath(input_gzi_index)
                )
            } else {
                gzindexes = gz_faidx(Channel.fromPath(margs.input_ref))
                gzindexes.ifEmpty{
                    if (params.containsKey("igv") && params.igv){
                        log.warn """\
                            The input reference is compressed but not with bgzip, which is required to create an index.
                            The workflow will proceed but it will not be possible to load the reference in the IGV Viewer.
                            To use the IGV Viewer, provide an uncompressed, or bgzip compressed version of the input reference next time you run the workflow.
                            """.stripIndent()
                    }
                }
            }
        } else {
            gzindexes = Channel.empty()
        }

        // Generate fai index if the file is either compressed, or if fai doesn't exists
        if (!is_compressed && file("${margs.input_ref}.fai").exists()){
            ref_idx = Channel.fromPath("${margs.input_ref}.fai")
        } else {
            ref_idx = faidx(ref)
        }

        // Generate CRAM cache
        if (margs.output_cache){
            cram_cache(ref)
            ref_cache = cram_cache.out.ref_cache
        } else {
            ref_cache = null
        }

        // Generate mmi index
        if (margs.output_mmi){
            ref_mmi = make_mmi(ref)
        } else {
            ref_mmi = null
        }

    emit:
        ref = ref
        ref_idx = ref_idx
        ref_cache = ref_cache
        ref_mmi = ref_mmi
        ref_gzidx = gzindexes
}
