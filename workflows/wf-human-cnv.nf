import groovy.json.JsonBuilder

include {
    getVersions;
    getParams;
    callCNV;
    makeReport
} from "../modules/local/wf-human-cnv.nf"

include { get_genome; } from '../modules/local/common.nf'

workflow cnv {
    take:
        bam_channel
        read_stats

    main:
        // truncate bam channel to remove meta to keep compat with cnv pipe
        bam = bam_channel.map{ it -> tuple(it[0], it[1]) }

        genome = get_genome(bam)
        genome_match_channel = genome.genome_build.ifEmpty{throw new Exception("The genome build detected in the BAM is not compatible with this workflow.)")}

        cnvs = callCNV(bam, genome.genome_build)

        software_versions = getVersions()
        workflow_params = getParams()

        report = makeReport(read_stats, cnvs, software_versions.collect(), workflow_params, genome.genome_build)

    emit:
        report
}