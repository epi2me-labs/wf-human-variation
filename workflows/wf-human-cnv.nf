import groovy.json.JsonBuilder

include {
    callCNV;
    getVersions;
    add_snp_tools_to_versions;
    getParams;
    bgzip_and_index_vcf;
    makeReport
} from "../modules/local/wf-human-cnv.nf"

include {
    mosdepth;
    annotate_vcf
} from "../modules/local/common.nf"

workflow cnv {
    take:
        bam
        ref
        clair3_vcf
        bed
    main:
        // get mosdepth results for window size 1000
        mosdepth(bam, bed, ref, "1000", false)
        mosdepth_stats = mosdepth.out.mosdepth_tuple
        mosdepth_summary = mosdepth.out.summary
        if (params.depth_intervals) {
            mosdepth_perbase = mosdepth.out.perbase
        } else {
            mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
        }

        mosdepth_all = mosdepth_stats.concat(mosdepth_summary).concat(mosdepth_perbase).collect()
        
        cnvs = callCNV(clair3_vcf, mosdepth_all, ref)
        spectre_vcf = cnvs.spectre_vcf
        spectre_vcf_bgzipped = bgzip_and_index_vcf(spectre_vcf)
        spectre_bed = cnvs.spectre_bed
        spectre_karyotype = cnvs.spectre_karyotype

        // check if SnpEff annotations have been requested
        if (!params.annotation) {
            spectre_final_vcf = spectre_vcf_bgzipped
        }
        else {
            // append '*' to indicate that annotation should be performed on all chr at once
            vcf_for_annotation = spectre_vcf_bgzipped.map{ it << '*' }
            spectre_final_vcf = annotate_vcf(vcf_for_annotation, "hg38", "cnv").annot_vcf
        }

        software_versions_tmp = getVersions()
        software_versions = add_snp_tools_to_versions(software_versions_tmp)
        workflow_params = getParams()
        if (params.output_report){
            report = makeReport(software_versions.collect(), workflow_params, spectre_bed, spectre_karyotype)
        } else {
            report = Channel.empty()
        }

    emit:
        spectre_final_vcf.concat(report)
}