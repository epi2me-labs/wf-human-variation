import groovy.json.JsonBuilder

process getParams {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "params.json"
    cache false
    cpus 1
    memory "2 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process configure_igv {
    publishDir "${params.out_dir}/", mode: 'copy', pattern: 'igv.json', enabled: params.containsKey("igv") && params.igv
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        // the python script will work out what to do with all the files based on their
        // extensions
        path "file-names.txt"
        val locus_str
        val aln_extra_opts
        val var_extra_opts
    output: path "igv.json"
    script:
    // the locus argument just makes sure that the initial view in IGV shows something
    // interesting
    String locus_arg = locus_str ? "--locus $locus_str" : ""
    // extra options for alignment tracks
    def aln_opts_json_str = \
        aln_extra_opts ? new JsonBuilder(aln_extra_opts).toPrettyString() : ""
    String aln_extra_opts_arg = \
        aln_extra_opts ? "--extra-alignment-opts extra-aln-opts.json" : ""
    // extra options for variant tracks
    def var_opts_json_str = \
        var_extra_opts ? new JsonBuilder(var_extra_opts).toPrettyString() : ""
    String var_extra_opts_arg = \
        var_extra_opts ? "--extra-vcf-opts extra-var-opts.json" : ""
    """
    # write out JSON files with extra options for the alignment and variant tracks
    echo '$aln_opts_json_str' > extra-aln-opts.json
    echo '$var_opts_json_str' > extra-var-opts.json

    workflow-glue configure_igv \
        --fofn file-names.txt \
        $locus_arg \
        $aln_extra_opts_arg \
        $var_extra_opts_arg \
    > igv.json
    """
}

