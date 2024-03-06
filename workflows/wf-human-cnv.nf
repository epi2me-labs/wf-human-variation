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
    mosdepth
} from "../modules/local/common.nf"

workflow cnv {
    take:
        bam
        ref
        clair3_vcf
        bed
        chromosome_codes
    main:
        // get mosdepth results for window size 1000
        mosdepth(bam, bed, ref, "1000")
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
        spectre_final_vcf = bgzip_and_index_vcf(spectre_vcf)
        spectre_bed = cnvs.spectre_bed

        software_versions_tmp = getVersions()
        software_versions = add_snp_tools_to_versions(software_versions_tmp)
        workflow_params = getParams()

        report = makeReport(software_versions.collect(), workflow_params, spectre_bed, chromosome_codes)

    emit:
        spectre_final_vcf.concat(report)
}