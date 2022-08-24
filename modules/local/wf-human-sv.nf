import groovy.json.JsonBuilder

process minimap2_ubam {
    label "wf_human_sv"
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        file reference
        file reads
    output:
        path "*.mm2.bam", emit: bam
        path "*.mm2.bam.bai", emit: bam_index
    script:
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${reads} | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${reference} - \
    | samtools sort -@ ${params.ubam_sort_threads} --write-index -o ${reads.baseName}.mm2.bam##idx##${reads.baseName}.mm2.bam.bai -
    """
}


process getAllChromosomesBed {
    label "wf_human_sv"
    cpus 1
    input:
        file reference
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    faidx --transform bed $reference > allChromosomes.bed
    """
}


// Remove unmapped (4), non-primary (256) and supplemental (2048) alignments
process filterBam {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
    output:
        path "*filtered.bam", emit: bam
        path "*filtered.bam.bai", emit: bam_index
    script:
        def name = bam.simpleName
    """
    samtools view -@ $task.cpus -F 2308 -o ${name}.filtered.bam ${bam}
    samtools index ${name}.filtered.bam
    """
}


// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
process sniffles2 {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
        file tr_bed
    output:
        path "*.sniffles.vcf", emit: vcf
    script:
        def name = bam.simpleName
        def tr_arg = tr_bed.name != 'OPTIONAL_FILE' ? "--tandem-repeats ${tr_bed}" : ''
    """
    sniffles \
        --threads $task.cpus \
        --sample-id $name \
        --output-rnames \
        --cluster-merge-pos $params.cluster_merge_pos \
        --input $bam \
        $tr_arg \
        --vcf ${name}.sniffles.vcf
    sed -i '/.:0:0:0:NULL/d' ${name}.sniffles.vcf
    """
}

process mosdepth {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
        file target_bed
    output:
        path "*.regions.bed.gz", emit: mosdepth_bed
        path "*.global.dist.txt", emit: mosdepth_dist
    script:
        def name = bam.simpleName
    """
    mosdepth \
        -x \
        -t $task.cpus \
        -b $target_bed \
        $name \
        $bam
    """
}


process filterCalls {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file mosdepth_bed
        file target_bed
    output:
        path "*.filtered.vcf", emit: vcf
    script:
        def name = vcf.simpleName
        def sv_types_joined = params.sv_types.split(',').join(" ")
    """
    get_filter_calls_command.py \
        --target_bedfile $target_bed \
        --vcf $vcf \
        --depth_bedfile $mosdepth_bed \
        --min_sv_length $params.min_sv_length \
        --max_sv_length $params.max_sv_length \
        --sv_types $sv_types_joined \
        --min_read_support $params.min_read_support \
        --min_read_support_limit $params.min_read_support_limit > filter.sh

    bash filter.sh > ${name}.filtered.vcf
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
        path "*.wf_sv.vcf", emit: vcf
    script:
        def name = vcf.simpleName
    """
    vcfsort $vcf > ${name}.wf_sv.vcf
    """
}


process indexVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "${vcf}.gz", emit: vcf_gz
        path "${vcf}.gz.tbi", emit: vcf_tbi
    """
    cat $vcf | bgziptabix ${vcf}.gz
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
    TRUVARI=\$(which truvari)
    python \$TRUVARI version | sed 's/ /,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
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
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file read_stats
        file read_depth
        file eval_json
        file versions
        path "params.json"
    output:
        path "wf-human-sv-*.html", emit: html
    script:
        def name = vcf.simpleName
        def report_name = "wf-human-sv-" + params.report_name + '.html'
        def readStats = read_stats ? "--reads_summary ${read_stats}" : ""
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
    """
    report_sv.py \
        $report_name \
        --vcf $vcf \
        $readStats \
        --read_depth $read_depth \
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
