import groovy.json.JsonBuilder

process callCNV {
    label "spectre"
    cpus 2
    memory 8.GB
    input:
        tuple path(vcf), path(vcf_index)
        path ("readstats/*")
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path ("spectre_output/${params.sample_name}.vcf"), emit: spectre_vcf
        path ("spectre_output/${params.sample_name}_cnv.bed"), emit: spectre_bed
        path ("spectre_output/expected_karyotype.json"), emit: spectre_karyotype
    script:
        """
        spectre CNVCaller \
        --bin-size 1000 \
        --threshhold-quantile 10 \
        --dist-proportion 0.3 \
        --coverage readstats/ \
        --sample-id ${params.sample_name} \
        --output-dir spectre_output/ \
        --reference ${ref} \
        --blacklist grch38_blacklist_0.3 \
        --min-cnv-len 80000 \
        --snv ${vcf} \
        --metadata grch38_metadata
        """
}

process bgzip_and_index_vcf {
    cpus 1
    input:
        path (spectre_vcf)
    output:
        tuple path("${params.sample_name}.wf_cnv.vcf.gz"), path("${params.sample_name}.wf_cnv.vcf.gz.tbi"), emit: spectre_final_vcf
    script:
        """
        bgzip ${spectre_vcf}
        mv ${spectre_vcf}.gz ${params.sample_name}.wf_cnv.vcf.gz
        tabix -f -p vcf ${params.sample_name}.wf_cnv.vcf.gz
        """
}

process getVersions {
    label "spectre"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        spectre version 2>&1 | awk '{print \$2","\$4}' >> versions.txt
        """
}

process add_snp_tools_to_versions {
    label "wf_human_snp"
    cpus 1
    memory "2 GB"
    input: path "old_versions.txt"
    output: path "versions.txt"
    script:
    """
    cp old_versions.txt versions.txt
    run_clair3.sh --version | sed 's/ /,/' >> versions.txt
    """
}

process getParams {
    label "wf_common"
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
        path "versions/*"
        path "params.json"
        path cnv_bed
        path karyotype_json
    output:
        path("*wf-human-cnv-report.html")
    script:
        def report_name = "${params.sample_name}.wf-human-cnv-report.html"
        """
        workflow-glue report_cnv_spectre \
            --sample_id ${params.sample_name} \
            --params params.json \
            --versions versions \
            --cnv_bed ${cnv_bed} \
            --karyotype_json ${karyotype_json} \
            -o $report_name \
            --workflow_version ${workflow.manifest.version}
        """

}
