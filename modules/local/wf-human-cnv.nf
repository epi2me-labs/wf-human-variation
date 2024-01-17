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
        path: { "${params.out_dir}/qdna_seq" },
        mode: 'copy',
        saveAs: { filename -> filename.toString() ==~ /.*vcf\.gz.*/ ? null : filename }
    ]
    input:
        tuple path(bam), path(bai)
        val(genome_build)
    output:
        tuple path("${params.sample_name}_combined.bed"), path("${params.sample_name}*"), path("${params.sample_name}_noise_plot.png"), path("${params.sample_name}_isobar_plot.png"), emit: cnv_output
        tuple path("${params.sample_name}.wf_cnv.vcf.gz"), path("${params.sample_name}.wf_cnv.vcf.gz.tbi"), emit: cnv_vcf
    script:
        """
        run_qdnaseq.r --bam ${bam} --out_prefix ${params.sample_name} --binsize ${params.bin_size} --reference ${genome_build}
        cut -f5 ${params.sample_name}_calls.bed | paste ${params.sample_name}_bins.bed - > ${params.sample_name}_combined.bed

        # VCF will be malformed if it contains one CNV (CW-1491), check and fix if necessary
        mv ${params.sample_name}_calls.vcf raw.vcf
        mv ${params.sample_name}_segs.vcf raw_segs.vcf
        fix_1491_vcf.py -i raw.vcf -o ${params.sample_name}.wf_cnv.vcf --sample_id ${params.sample_name}
        fix_1491_vcf.py -i raw_segs.vcf -o ${params.sample_name}_segs.vcf --sample_id ${params.sample_name}

        # bgzip and index calls VCF
        bgzip ${params.sample_name}.wf_cnv.vcf
        tabix -f -p vcf ${params.sample_name}.wf_cnv.vcf.gz
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
    cpus 1
    memory 12.GB
    input:
        path(read_stats)
        tuple path(cnv_calls), val(cnv_files), path(noise_plot), path(isobar_plot)
        path "versions/*"
        path "params.json"
        val(genome_build)
    output:
        path("*wf-human-cnv-report.html")

    script:
        def report_name = "${params.sample_name}.wf-human-cnv-report.html"
        """
        workflow-glue cnv_plot \
            -q ${cnv_calls} \
            -o $report_name \
            --read_stats ${read_stats}\
            --params params.json \
            --versions versions \
            --bin_size ${params.bin_size} \
            --genome ${genome_build} \
            --sample_id ${params.sample_name} \
            --noise_plot ${noise_plot} \
            --isobar_plot ${isobar_plot}
        """

}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_cnv {
    // publish inputs to output directory
    label "wf_cnv"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


