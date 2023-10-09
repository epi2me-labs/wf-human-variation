import groovy.json.JsonBuilder

process call_str {
    // first subset the repeats BED file, then use this to do straglr genotyping
    label "wf_human_str"
    cpus 1
    input:
        tuple val(chr), path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path(repeat_bed)
    output:
        tuple path("*.vcf.gz"), path ("*.tsv"), val(chr), optional: true, emit: straglr_output
    shell:
        '''
        { grep !{chr} -Fw !{repeat_bed} || true; } > repeats_subset.bed

        if [[ -s repeats_subset.bed ]]; then
            straglr-genotype --loci repeats_subset.bed \
                --sample !{params.sample_name} \
                --tsv !{chr}_straglr.tsv \
                -v !{chr}_tmp.vcf \
                --sex !{params.sex} !{xam} !{ref} \
                --min_support 1 \
                --min_cluster_size 1
            cat !{chr}_tmp.vcf | vcfstreamsort | bgziptabix !{chr}_straglr.vcf.gz
        else
            echo "blank subset BED"
        fi
        '''
}


process annotate_repeat_expansions {
    // annotate using Stranger
    label "wf_human_str"
    input:
        tuple path(vcf), path (tsv), val (chr)
        path(variant_catalogue_hg38)
    output:
        tuple path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: annotated_vcf
        path "*.tsv", emit: annotated_tsv
        path "*_plot.tsv", emit: plot_tsv
        path "*_annotated.tsv", emit: stranger_annotation
    script:
        """
        stranger -f ${variant_catalogue_hg38} ${vcf} \
            | sed 's/\\ /_/g' \
            | bgziptabix ${chr}_repeat-expansion_annotated.vcf.gz
        SnpSift extractFields ${chr}_repeat-expansion_annotated.vcf.gz \
            CHROM POS ALT FILTER REF RL RU REPID VARID STR_STATUS > ${chr}_repeat-expansion_annotated.tsv
        SnpSift extractFields ${chr}_repeat-expansion_annotated.vcf.gz \
            CHROM POS DisplayRU STR_NORMAL_MAX STR_PATHOLOGIC_MIN VARID Disease > ${chr}_repeat-expansion_plot.tsv
        """
}


process merge_tsv {
    // merge the contig TSVs
    input:
        path (plot_tsv)
        path (straglr_tsv)
        path (stranger_tsv)
    output:
        tuple path ("*plot.tsv"), path ("*straglr.tsv"), path ("*stranger.tsv")
    script:
        """
        awk 'NR == 1 || FNR > 1' ${plot_tsv} >${params.sample_name}_plot.tsv
        awk 'NR == 1 || FNR > 1' ${stranger_tsv} >${params.sample_name}_stranger.tsv
        # remove the first line from straglr_tsv as it has a two-line header
        sed -i '1d' ${straglr_tsv}
        awk 'NR == 1 || FNR > 1' ${straglr_tsv} > ${params.sample_name}.wf_str.straglr.tsv
        """
}


process merge_vcf {
    // merge the contig VCFs
    input:
        path (vcfs)
        path (vcf_indexes)
    output:
        tuple path ("*str.vcf.gz"), path("*str.vcf.gz.tbi"), emit: final_vcf
    script:
        """
        bcftools concat ${vcfs} > ${params.sample_name}.wf_str.vcf
        bgzip ${params.sample_name}.wf_str.vcf
        tabix ${params.sample_name}.wf_str.vcf.gz
        """
}


process make_report {
    input:
        path(vcf)
        path(straglr_tsv) 
        path(plot_tsv)
        path(stranger_annotation)
        path "versions.txt"
        path "params.json"
        path(bam_stats)
    output:
        path "*wf-human-str-report.html", emit: html
    script:
        def report_name = "${params.sample_name}.wf-human-str-report.html"
        """
        workflow-glue report_str \
            -o $report_name \
            --params params.json \
            --sample_name ${params.sample_name} \
            --version versions.txt \
            --vcf ${vcf} \
            --straglr ${straglr_tsv} \
            --stranger ${plot_tsv} \
            --stranger_annotation ${stranger_annotation} \
            --read_stats ${bam_stats}
        """
}


process getVersions {
    label "wf_human_str"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
        samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
        straglr-genotype --version | sed 's/^/straglr-genotype,/' >> versions.txt
        stranger --version | sed 's/^/stranger,/' >> versions.txt
        """
}


process getParams {
    label "wf_human_str"
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


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_str {
    // publish inputs to output directory
    label "wf_human_str"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}
