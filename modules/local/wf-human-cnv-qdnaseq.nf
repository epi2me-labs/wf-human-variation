import groovy.json.JsonBuilder

// fix_vcf was unglued to avoid installing base deps in CNV container
process callCNV {
    label "wf_cnv"
    cpus 1
    memory { 16.GB * task.attempt }
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    // publish everything except the cnv_vcf to qdna_seq directory
    publishDir = [
        path: { "${params.outdirNonPhi}/qdna_seq" },
        mode: 'copy',
        saveAs: { filename -> filename.toString() ==~ /.*vcf\.gz.*/ ? null : filename }
    ]
    input:
        tuple path(bam), path(bai), val(xam_meta)
        val(genome_build)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}_combined.bed"), path("${xam_meta.alias}*"), path("${xam_meta.alias}_noise_plot.png"), path("${xam_meta.alias}_isobar_plot.png"), emit: cnv_output
        tuple val(xam_meta), path("${xam_meta.alias}.wf_cnv.vcf.gz"), path("${xam_meta.alias}.wf_cnv.vcf.gz.tbi"), emit: cnv_vcf
    script:
        """
        run_qdnaseq.r --bam ${bam} --out_prefix ${xam_meta.alias} --binsize ${params.qdnaseq_bin_size} --reference ${genome_build}
        cut -f5 ${xam_meta.alias}_calls.bed | paste ${xam_meta.alias}_bins.bed - > ${xam_meta.alias}_combined.bed

        # VCF will be malformed if it contains one CNV (CW-1491), check and fix if necessary
        mv ${xam_meta.alias}_calls.vcf raw.vcf
        mv ${xam_meta.alias}_segs.vcf raw_segs.vcf
        fix_1491_vcf.py -i raw.vcf -o ${xam_meta.alias}.wf_cnv.vcf --sample_id ${xam_meta.alias}
        fix_1491_vcf.py -i raw_segs.vcf -o ${xam_meta.alias}_segs.vcf --sample_id ${xam_meta.alias}

        # bgzip and index calls VCF
        bgzip ${xam_meta.alias}.wf_cnv.vcf
        tabix -f -p vcf ${xam_meta.alias}.wf_cnv.vcf.gz
        """
}

process getVersions {
    label "wf_cnv"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
        R --version | grep -w R | grep version | cut -f3 -d" " | sed 's/^/R,/' >> versions.txt
        R --slave -e 'packageVersion("QDNAseq")' | cut -d\\' -f2 | sed 's/^/QDNAseq,/' >> versions.txt
        samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
        """
}


process getParams {
    label "wf_cnv"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        """
        # Output nextflow params object to JSON
        echo '$paramsJSON' > params.json
        """
}


process makeReport {
    label "wf_common"
    cpus 1
    memory 12.GB
    input:
        path(read_stats)
        tuple val(xam_meta), path(cnv_calls), val(cnv_files), path(noise_plot), path(isobar_plot)
        path "versions/*"
        path "params.json"
        val(genome_build)
    output:
        path("*wf-human-cnv-report.html")

    script:
        def report_name = "${xam_meta.alias}.wf-human-cnv-report.html"
        """
        workflow-glue report_cnv_qdnaseq \
            -q ${cnv_calls} \
            -o $report_name \
            --read_stats ${read_stats}\
            --params params.json \
            --versions versions \
            --bin_size ${params.qdnaseq_bin_size} \
            --genome ${genome_build} \
            --sample_id ${xam_meta.alias} \
            --noise_plot ${noise_plot} \
            --isobar_plot ${isobar_plot} \
            --workflow_version ${workflow.manifest.version}
        """

}


