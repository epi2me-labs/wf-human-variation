import groovy.json.JsonBuilder

process getParams {
    label "wf_common"
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
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        // the python script will work out what to do with all the files based on their
        // extensions
        path "file-names.txt"
        val locus_str
        val bam_extra_opts
        val vcf_extra_opts
    output: path "igv.json"
    script:
    // the locus argument just makes sure that the initial view in IGV shows something
    // interesting
    String locus_arg = locus_str ? "--locus $locus_str" : ""
    // extra options for alignment tracks
    def bam_opts_json_str = \
        bam_extra_opts ? new JsonBuilder(bam_extra_opts).toPrettyString() : ""
    String bam_extra_opts_arg = \
        bam_extra_opts ? "--extra-bam-opts bam-extra-opts.json" : ""
    // extra options for variant tracks
    def vcf_opts_json_str = \
        vcf_extra_opts ? new JsonBuilder(vcf_extra_opts).toPrettyString() : ""
    String vcf_extra_opts_arg = \
        vcf_extra_opts ? "--extra-vcf-opts vcf-extra-opts.json" : ""
    """
    # write out JSON files with extra options for the alignment and variant tracks
    echo '$bam_opts_json_str' > bam-extra-opts.json
    echo '$vcf_opts_json_str' > vcf-extra-opts.json

    workflow-glue configure_igv \
        --fofn file-names.txt \
        $locus_arg \
        $bam_extra_opts_arg \
        $vcf_extra_opts_arg \
    > igv.json
    """
}

