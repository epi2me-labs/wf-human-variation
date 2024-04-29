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
        // truncate bam channel to remove meta to keep compat with cnv pipe
        bam = bam_channel.map{ it -> tuple(it[0], it[1]) }

        cnvs = callCNV(bam, genome_build)

        software_versions = getVersions()
        workflow_params = getParams()


        if (params.output_report){
            report = makeReport(read_stats, cnvs.cnv_output, software_versions.collect(), workflow_params, genome_build)
        } else {
            report = Channel.empty()
        }

    emit:
        cnvs.cnv_vcf.concat(report)
}