include {
    concat_vcfs as concat_str_vcfs;
} from "../modules/local/common.nf"
include {
    call_str;
    annotate_repeat_expansions;
    merge_tsv;
    bam_region_filter;
    bam_read_filter;
    generate_str_content;
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

    // subset contig BAM to include only STR regions
    // ignore those contigs which aren't in repeats BED so this output is optional
    str_regions_bam = bam_region_filter(bam_channel, str_list)

    // make a channel of chr, xam, xam_idx, vcf, straglr_tsv
    reads_join = str_regions_bam.join(str_vcf_and_tsv.straglr_output)

    // subset contig STR regions BAM to include only supporting reads from straglr 
    str_reads_bam = bam_read_filter(reads_join)

    // join all the contig information ready to generate STR content
    str_content_join = str_vcf_and_tsv.straglr_output.join(annotations.stranger_annotation).join(str_reads_bam.reads_bam)
    
    str_content = generate_str_content(
      str_content_join,
      str_list
    )

    // merge the contig VCFs
    stranger_vcfs_and_tbis = annotations.stranger_annotation.map{ it -> [it[1], it[2]] }.collect()
    merged_vcf = concat_str_vcfs(stranger_vcfs_and_tbis, "${params.sample_name}.wf_str").final_vcf

    // merge the contig TSVs/CSVs
    plot_tsv_all = annotations.stranger_annotation.map{it -> it[3]}.collect()
    straglr_tsv_all = str_vcf_and_tsv.map{it -> it[2]}.collect()
    stranger_annotation_all = annotations.stranger_annotation.map{it -> it[4]}.collect()

    str_content_all = str_content.str_content.collect()

    // get the merged TSVs ready for the report
    merged_tsvs = merge_tsv(plot_tsv_all, straglr_tsv_all, stranger_annotation_all, str_content_all)
    merged_plot = merged_tsvs.map{it -> it[0]}
    merged_straglr = merged_tsvs.map{it -> it[1]}
    merged_stranger = merged_tsvs.map{it -> it[2]}
    merged_str_content = merged_tsvs.map{it -> it[3]}

    report = make_report(merged_vcf, merged_straglr, merged_plot, merged_stranger, merged_str_content, software_versions, workflow_params, read_stats)

  emit:
    merged_vcf.concat(report).concat(merged_straglr).flatten()

}
