// This subworkflow aims at pre-processing an input reference
// genome, preparing the appropriate indexes and cache.
// By default, the workflow takes an input reference,
// decompress it if it is gzipped, and import (if available)
// or index with `samtools faidx` the reference genome.
// if the user requests it, the workflow can also:
// 1. generate the cram cache, with annex environmental
//    `REF_PATH` variable
// 2. generate the minimap2 `.mmi` index, for faster alignment

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
            ref = decompress_ref(ref)
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
}
