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
    sex

  main:
    // turn ref channel into value channel so it can be used more than once
    ref_as_value = ref_channel.collect()

    str_list = projectDir.resolve("./data/wf_str_repeats.bed").toString()
    variant_catalogue_hg38 = projectDir.resolve("./data/variant_catalog_hg38.json").toString()

    // call straglr and get annotations per contig
    str_vcf_and_tsv = call_str(bam_channel.combine(sex), ref_as_value, str_list)
    annotations = annotate_repeat_expansions(str_vcf_and_tsv, variant_catalogue_hg38)

    software_versions = getVersions()
    workflow_params = getParams()

    // subset contig BAM to include only STR regions
    // ignore those contigs which aren't in repeats BED so this output is optional
    str_regions_bam = bam_region_filter(bam_channel, str_list) | map { bam, bai, meta -> tuple(meta.sq, bam, bai, meta) }

    // make a channel of chr, xam, xam_idx, vcf, straglr_tsv
    reads_join = str_regions_bam.join(str_vcf_and_tsv)

    // subset contig STR regions BAM to include only supporting reads from straglr
    str_reads_bam = bam_read_filter(reads_join) | map { bam, bai, meta -> tuple(meta.sq, bam, bai, meta) }

    // join all the contig information ready to generate STR content
    str_content_join = str_vcf_and_tsv.join(annotations).join(str_reads_bam)

    str_content = generate_str_content(
      str_content_join,
      str_list
    ).transpose().groupTuple()


    branched_annotations = annotations.multiMap { chr, vcf, tbi, plot, annot ->
        stranger_vcfs_and_tbis: vcf
        plot_tsv_all: plot
        stranger_annotations: annot
    }

    // merge the contig VCFs
    // bam_channel.xam_meta is per contig (containing sq: and id:) 
    // so just use meta.alias to combine with the grouped branched_annotation vcfs
    merged_vcf = concat_str_vcfs(
        bam_channel.map{xam, xai, meta ->   ['alias': meta.alias] }.unique().
        combine(branched_annotations.stranger_vcfs_and_tbis).groupTuple(),
        "wf_str"
    ).final_vcf

    // merge the contig TSVs/CSVs
    straglr_tsv_all = str_vcf_and_tsv.map{ chr, vcf, tsv -> tsv }.collect()
    branched_merged = merge_tsv(
        branched_annotations.plot_tsv_all.collect(),
        straglr_tsv_all.collect(),
        branched_annotations.stranger_annotations.collect(),
        str_content)
        | multiMap { plot, straglr_table, stranger_table, str_content_table ->
            plot: plot
            straglr: straglr_table
            stranger: stranger_table
            str_content: str_content_table
        }

    if (params.output_report){
      report = make_report(
          merged_vcf,
          branched_merged.straglr.collect(),
          branched_merged.plot.collect(),
          branched_merged.stranger.collect(),
          branched_merged.str_content.collect(),
          software_versions,
          workflow_params,
          read_stats,
          sex
      )
    } else {
      report = Channel.empty()
    }

  emit:
    output = merged_vcf.map{meta, vcf, tbi -> [vcf, tbi]}.concat(report).concat(branched_merged.straglr).flatten()
    str_vcf = merged_vcf.collect()

}
