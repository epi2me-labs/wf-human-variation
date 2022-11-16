import groovy.json.JsonBuilder

include {
    filterBam;
    sniffles2;
    filterCalls;
    sortVCF;
    indexVCF;
    getVersions;
    getParams;
    report;
} from "../modules/local/wf-human-sv.nf"


workflow bam {
    take:
        bam_channel
        reference
        target
        mosdepth_stats
        optional_file
    main:
        // truncate bam channel to remove meta to keep compat with sv pipe
        bam = bam_channel.map{ it -> tuple(it[0], it[1]) }

        called = variantCall(bam, reference, target, mosdepth_stats, optional_file)
        report = runReport(
            called.vcf.collect(),
            [],
            mosdepth_stats,
            optional_file)
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
        )
}


workflow variantCall {
    take:
        bam
        reference
        target_bed
        mosdepth_stats
        optional_file
    main:

        // tandom_repeat bed
        if(params.tr_bed == null) {
            tr_bed = optional_file
        } else {
            tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true)
        }

        filterBam(bam, reference)
        sniffles2(filterBam.out.cram, tr_bed, reference)
        filterCalls(sniffles2.out.vcf, mosdepth_stats, target_bed)
        sortVCF(filterCalls.out.vcf)
        indexVCF(sortVCF.out.vcf)
    emit:
        vcf = indexVCF.out.vcf_gz
        vcf_index = indexVCF.out.vcf_tbi
}

workflow runReport {
    take:
        vcf
        read_stats
        depth_bed
        eval_json
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        report(
            vcf.collect(),
            read_stats.collect(),
            depth_bed.collect(),
            eval_json,
            software_versions, 
            workflow_params)
    emit:
        html = report.out.html
}

