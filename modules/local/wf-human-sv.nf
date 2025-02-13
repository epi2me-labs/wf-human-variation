import groovy.json.JsonBuilder

// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
// --input-exclude-flags 2308: Remove unmapped (4), non-primary (256) and supplemental (2048) alignments
process sniffles2 {
    label "wf_human_sv"
    cpus params.threads
    memory 6.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        file tr_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val genome_build
    output:
        tuple val(xam_meta), path("*.sniffles.vcf"), emit: vcf
        path "${xam_meta.alias}.wf_sv.snf", emit: snf
    publishDir \
        path: "${params.outdirNonPhi}",
        pattern: "${xam_meta.alias}.wf_sv.snf",
        mode: 'copy'
    script:
        // if tr_arg is not provided and genome_build is set
        // automatically pick the relevant TR BED from the SV image
        def tr_arg = ""
        if (tr_bed.name != 'OPTIONAL_FILE'){
            tr_arg = "--tandem-repeats ${tr_bed}"
        }
        else if (genome_build) {
            log.warn "Automatically selecting TR BED: ${genome_build}.trf.bed"
            tr_arg = "--tandem-repeats \${WFSV_TRBED_PATH}/${genome_build}.trf.bed"
        }
        def sniffles_args = params.sniffles_args ?: ''
        def min_sv_len = params.min_sv_length ? "--minsvlen ${params.min_sv_length}" : ""
        // Perform internal phasing only if snp not requested; otherwise, use joint phasing.
        def phase = params.phased ? "--phase" : ""
    """
    sniffles \
        --threads $task.cpus \
        --sample-id ${xam_meta.alias} \
        --output-rnames \
        ${min_sv_len} \
        --cluster-merge-pos $params.cluster_merge_pos \
        --input $xam \
        --reference $ref \
        --input-exclude-flags 2308 \
        --snf ${xam_meta.alias}.wf_sv.snf \
        $tr_arg \
        $sniffles_args \
        $phase \
        --vcf ${xam_meta.alias}.sniffles.vcf
    sed '/.:0:0:0:NULL/d' ${xam_meta.alias}.sniffles.vcf > tmp.vcf
    mv tmp.vcf ${xam_meta.alias}.sniffles.vcf
    """
}


process filterCalls {
    cpus { params.threads < 2 ? 2 : params.threads }
    memory 4.GB
    input:
        tuple val(xam_meta), path(vcf)
        path mosdepth_summary // MOSDEPTH_TUPLE
        file target_bed
        val chromosome_codes
    output:
        tuple val(xam_meta), path("*.filtered.vcf"), emit: vcf
    script:
    String ctgs = chromosome_codes.join(',')
    def ctgs_filter = params.include_all_ctgs ? "" : "--contigs ${ctgs}"
    """
    # Filter contigs requre the input VCF to be compressed and indexed
    bcftools view -O z $vcf > input.vcf.gz && tabix -p vcf input.vcf.gz

    # Create filtering script
    get_filter_calls_command.py \
        --bcftools_threads $task.cpus \
        --target_bedfile $target_bed \
        --vcf input.vcf.gz \
        --depth_summary $mosdepth_summary \
        --min_read_support $params.min_read_support \
        --min_read_support_limit $params.min_read_support_limit \
        ${ctgs_filter} > filter.sh

    # Run filtering
    bash filter.sh > ${xam_meta.alias}.filtered.vcf
    """
}


// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
process sortVCF {
    label "wf_human_sv"
    cpus 2
    memory 4.GB
    input:
        tuple val(xam_meta), path(vcf)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_sv.vcf.gz"), emit: vcf_gz
        tuple val(xam_meta), path("${xam_meta.alias}.wf_sv.vcf.gz.tbi"), emit: vcf_tbi
    script:
    """
    bcftools sort -m 2G -T ./ -O z $vcf > ${xam_meta.alias}.wf_sv.vcf.gz
    tabix -p vcf ${xam_meta.alias}.wf_sv.vcf.gz
    """
}


process getVersions {
    label "wf_human_sv"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    trap '' PIPE # suppress SIGPIPE without interfering with pipefail
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    truvari version | sed 's/ /,/' >> versions.txt
    sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    echo `seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d ':' -f 2 | sed 's/ /seqtk,/'` >> versions.txt
    """
}


process getParams {
    label "wf_human_sv"
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


process report {
    label "wf_common"
    cpus 1
    memory 6.GB
    input:
        tuple val(xam_meta), path(vcf)
        file eval_json
        file versions
        path "params.json"
    output:
        path "*report.html", emit: html, optional: true
        path "${xam_meta.alias}.svs.json", emit: json
    script:
        def report_name = "${xam_meta.alias}.wf-human-sv-report.html"
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
        def generate_html = params.output_report ? "" : "--skip_report"
    """
    workflow-glue report_sv \
        $report_name \
        --vcf $vcf \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        --output_json "${xam_meta.alias}.svs.json" \
        --workflow_version ${workflow.manifest.version} \
        $evalResults $generate_html
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_sv {
    // publish inputs to output directory
    label "wf_human_sv"
    publishDir "${params.outdirNonPhi}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}
