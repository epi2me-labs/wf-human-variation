import groovy.json.JsonBuilder

include {
    getVersions;
    getParams;
    callCNV;
    makeReport
} from "../modules/local/wf-human-cnv-qdnaseq.nf"


workflow cnv {
    take:
        bam_channel
        read_stats
        genome_build

    main:

        cnvs = callCNV(bam_channel, genome_build)

        software_versions = getVersions()
        workflow_params = getParams()


        if (params.output_report){
            report = makeReport(read_stats, cnvs.cnv_output, software_versions.collect(), workflow_params, genome_build)
        } else {
            report = Channel.empty()
        }

    emit:
        output = cnvs.cnv_vcf.map{ meta, vcf, tbi -> [vcf, tbi] }.concat(report)
        cnv_vcf = cnvs.cnv_vcf
}
