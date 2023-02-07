import groovy.json.JsonBuilder

include {
    getGenome;
    getVersions;
    getParams;
    callCNV;
    makeReport
} from "../modules/local/wf-human-cnv.nf"


workflow cnv {
    take:
        bam_channel
        read_stats

    main:
        // truncate bam channel to remove meta to keep compat with cnv pipe
        bam = bam_channel.map{ it -> tuple(it[0], it[1]) }

        genome = getGenome(bam)
        genome_match_channel = genome.genome_build.ifEmpty{exit 1, log.error('The genome build detected in the BAM is not compatible with this workflow.')}

        cnvs = callCNV(bam, genome.genome_build)

        software_versions = getVersions()
        workflow_params = getParams()

        report = makeReport(read_stats, cnvs, software_versions.collect(), workflow_params, genome.genome_build)

    emit:
        report
}