import groovy.json.JsonBuilder

process callCNV {
    label "spectre"
    cpus 2
    memory 8.GB
    input:
        tuple val(xam_meta), path(vcf), path(vcf_index)
        path("readstats/*")
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val(genome_build)
    output:
        tuple val(xam_meta), path("spectre_output/${xam_meta.alias}.vcf"), emit: spectre_vcf
        tuple val(xam_meta), path("spectre_output/${xam_meta.alias}_cnv.bed"), emit: spectre_bed
        tuple val(xam_meta), path("spectre_output/predicted_karyotype.txt"), emit: spectre_karyotype
    script:
        def spectre_args = params.spectre_args ?: ''
        """
        spectre CNVCaller \
        --bin-size 1000 \
        --coverage readstats/ \
        --sample-id ${xam_meta.alias} \
        --output-dir spectre_output/ \
        --reference ${ref} \
        --blacklist ${genome_build}_blacklist_v1.0 \
        --snv ${vcf} \
        --metadata ${genome_build}_metadata \
        $spectre_args
        """
}

process bgzip_and_index_vcf {
    cpus 1
    input:
        tuple val(xam_meta), path(spectre_vcf)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_cnv.vcf.gz"), path("${xam_meta.alias}.wf_cnv.vcf.gz.tbi"), emit: spectre_final_vcf
    script:
        """
        bgzip ${spectre_vcf}
        mv ${spectre_vcf}.gz ${xam_meta.alias}.wf_cnv.vcf.gz
        tabix -f -p vcf ${xam_meta.alias}.wf_cnv.vcf.gz
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


process makeReport {
    label "wf_common"
    cpus 1
    memory 12.GB
    input:
        path "versions/*"
        path "params.json"
        tuple val(xam_meta), path(cnv_bed)
        tuple val(xam_meta), path(karyotype)
        val genome_build
    output:
        path("*wf-human-cnv-report.html")
    script:
        def report_name = "${xam_meta.alias}.wf-human-cnv-report.html"
        """
        workflow-glue report_cnv_spectre \
            --sample_id ${xam_meta.alias} \
            --params params.json \
            --versions versions \
            --cnv_bed ${cnv_bed} \
            --karyotype ${karyotype} \
            -o $report_name \
            --workflow_version ${workflow.manifest.version} \
            --genome_build $genome_build
        """

}
