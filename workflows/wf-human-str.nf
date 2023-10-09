include {
    call_str;
    annotate_repeat_expansions;
    merge_vcf;
    merge_tsv;
    make_report;
    getVersions;
    getParams
} from "../modules/local/wf-human-str.nf"

// workflow module
workflow str {
  take:
    bam_channel
    ref_channel
    read_stats

  main:
    // turn ref channel into value channel so it can be used more than once
    ref_as_value = ref_channel.collect()

    str_list = projectDir.resolve("./data/wf_str_repeats.bed").toString()
    variant_catalogue_hg38 = projectDir.resolve("./data/variant_catalog_hg38.json").toString()

    // call straglr and get annotations per contig
    str_vcf_and_tsv = call_str(bam_channel, ref_as_value, str_list)

    annotations = annotate_repeat_expansions(str_vcf_and_tsv.straglr_output, variant_catalogue_hg38)

    software_versions = getVersions()
    workflow_params = getParams()

    // merge the contig VCFs
    annotations.annotated_vcf.collect{it[0]}.set { all_vcfs }
    annotations.annotated_vcf.collect{it[1]}.set { all_vcf_indexes }
    merged_vcf_and_index = merge_vcf(all_vcfs, all_vcf_indexes)
    merged_vcf_and_index.collect{it[0]}.set { merged_vcf }
  
    // merge the contig TSVs
    plot_tsv_all = annotations.plot_tsv.collect()
    straglr_tsv_all = str_vcf_and_tsv.map{it -> it[1]}.collect()
    stranger_annotation_all = annotations.stranger_annotation.collect()

    // get the merged TSVs ready for the report
    merged_tsvs = merge_tsv(plot_tsv_all, straglr_tsv_all, stranger_annotation_all)
    merged_plot = merged_tsvs.map{it -> it[0]}
    merged_straglr = merged_tsvs.map{it -> it[1]}
    merged_stranger = merged_tsvs.map{it -> it[2]}

    report = make_report(merged_vcf, merged_straglr, merged_plot, merged_stranger, software_versions, workflow_params, read_stats)

  emit:
    merged_vcf_and_index.concat(report).concat(merged_straglr).flatten()

}
