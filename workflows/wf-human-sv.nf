import groovy.json.JsonBuilder

include {
    getAllChromosomesBed;
    filterBam;
    sniffles2;
    mosdepth;
    filterCalls;
    sortVCF;
    indexVCF;
    getVersions;
    getParams;
    report;
} from "../modules/local/wf-human-sv.nf"

workflow bam {
    take:
        bam
        bam_index
        reference
        target
        optional_file
    main:
        if(!target) {
            target = getAllChromosomesBed(reference).all_chromosomes_bed
        }

        called = variantCall(bam, bam_index, reference, target, optional_file)
        report = runReport(
            called.vcf.collect(),
            [],
            called.read_depth.collect(),
            optional_file)
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
            bam,
            bam_index,
        )
}


workflow variantCall {
    take:
        bam
        bai
        reference
        target_bed
        optional_file
    main:

        // tandom_repeat bed
        if(params.tr_bed == null) {
            tr_bed = optional_file
        } else {
            tr_bed = Channel.fromPath(params.tr_bed, checkIfExists: true)
        }

        filterBam(bam, bai)
        sniffles2(filterBam.out.bam, filterBam.out.bam_index, tr_bed)
        mosdepth(filterBam.out.bam, filterBam.out.bam_index, target_bed)
        filterCalls(sniffles2.out.vcf, mosdepth.out.mosdepth_bed, target_bed)
        sortVCF(filterCalls.out.vcf)
        indexVCF(sortVCF.out.vcf)
    emit:
        vcf = indexVCF.out.vcf_gz
        vcf_index = indexVCF.out.vcf_tbi
        read_depth = mosdepth.out.mosdepth_dist
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

