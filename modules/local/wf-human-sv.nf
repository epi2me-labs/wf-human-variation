import groovy.json.JsonBuilder

// Remove unmapped (4), non-primary (256) and supplemental (2048) alignments
process filterBam {
    label "wf_human_sv"
    cpus params.threads
    input:
        tuple path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(xam_fmt), val(xai_fmt)
    output:
        tuple path("${params.sample_name}.filtered.${xam_fmt}"), path("${params.sample_name}.filtered.${xam_fmt}.${xai_fmt}"), emit: xam
    script:
    """
    samtools view -@ $task.cpus -F 2308 -o ${params.sample_name}.filtered.${xam_fmt}##idx##${params.sample_name}.filtered.${xam_fmt}.${xai_fmt} -O ${xam_fmt} --write-index --reference ${ref} ${xam}
    """
}


// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
process sniffles2 {
    label "wf_human_sv"
    cpus params.threads
    input:
        tuple path(xam), path(xam_idx)
        file tr_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "*.sniffles.vcf", emit: vcf
    script:
        def tr_arg = tr_bed.name != 'OPTIONAL_FILE' ? "--tandem-repeats ${tr_bed}" : ''
        def sniffles_args = params.sniffles_args ?: ''
        def min_sv_len = params.min_sv_length ? "--minsvlen ${params.min_sv_length}" : ""
        // Perform internal phasing only if snp not requested; otherwise, use joint phasing.
        def phase = params.phased && (!params.snp || params.output_separate_phased) ? "--phase" : ""
    """
    sniffles \
        --threads $task.cpus \
        --sample-id ${params.sample_name} \
        --output-rnames \
        ${min_sv_len} \
        --cluster-merge-pos $params.cluster_merge_pos \
        --input $xam \
        $tr_arg \
        $sniffles_args \
        $phase \
        --vcf ${params.sample_name}.sniffles.vcf
    sed '/.:0:0:0:NULL/d' ${params.sample_name}.sniffles.vcf > tmp.vcf
    mv tmp.vcf ${params.sample_name}.sniffles.vcf
    """
}


process filterCalls {
    cpus 1
    input:
        file vcf
        path mosdepth_summary // MOSDEPTH_TUPLE
        file target_bed
    output:
        path "*.filtered.vcf", emit: vcf
    script:
    """
    get_filter_calls_command.py \
        --target_bedfile $target_bed \
        --vcf $vcf \
        --depth_summary $mosdepth_summary \
        --min_read_support $params.min_read_support \
        --min_read_support_limit $params.min_read_support_limit > filter.sh

    bash filter.sh > ${params.sample_name}.filtered.vcf
    """
}


// NOTE This is the last touch the VCF has as part of the workflow,
//  we'll rename it with its desired output name here
process sortVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "${params.sample_name}.wf_sv.vcf.gz", emit: vcf_gz
        path "${params.sample_name}.wf_sv.vcf.gz.tbi", emit: vcf_tbi
    script:
    """
    bcftools sort -m 2G -T ./ -O z $vcf > ${params.sample_name}.wf_sv.vcf.gz
    tabix -p vcf ${params.sample_name}.wf_sv.vcf.gz
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
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
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
    cpus 1
    input:
        file vcf
        file eval_json
        file versions
        path "params.json"
    output:
        path "*report.html", emit: html
    script:
        def report_name = "${params.sample_name}.wf-human-sv-report.html"
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
    """
    workflow-glue report_sv \
        $report_name \
        --vcf $vcf \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        $evalResults
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_sv {
    // publish inputs to output directory
    label "wf_human_sv"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}
