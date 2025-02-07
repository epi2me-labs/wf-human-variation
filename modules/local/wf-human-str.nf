import groovy.json.JsonBuilder

process call_str {
    // first subset the repeats BED file, then use this to do straglr genotyping
    // sex comes from either params.sex or from inferred_sex
    //   and is assumed to be a non-null value of either XX or XY
    //   which will need translating to female or male for straglr
    label "wf_human_str"
    // Occasionally, straglr shows 150% of usage; give two cores to prevent the issue. 
    cpus 2
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta), val(sex)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path(repeat_bed)
    output:
        // emit meta.sq as join key (TODO a multi sample approach will want a compound key)
        tuple val(xam_meta.sq), path("*.vcf.gz"), path ("*.tsv"), optional: true
    script:
        // workflow default behaviour is to assume XX, so fall back if sex is not XY
        String straglr_sex = sex == "XY" ? "male" : "female"
        def chr = xam_meta.sq
        """
        { grep '${chr}' -Fw ${repeat_bed} || true; } > repeats_subset.bed
        if [[ -s repeats_subset.bed ]]; then
            straglr-genotype --loci repeats_subset.bed \
                --sample ${xam_meta.alias} \
                --tsv ${chr}_straglr.tsv \
                -v ${chr}_tmp.vcf \
                --sex ${straglr_sex} \
                --min_support 1 \
                --threads 1 \
                --min_cluster_size 1 \
                ${xam} ${ref}
            bgzip -c ${chr}_tmp.vcf > ${chr}_straglr.vcf.gz
            tabix -p vcf ${chr}_straglr.vcf.gz
        else
            echo "blank subset BED"
        fi
        """
}


process annotate_repeat_expansions {
    // annotate using Stranger
    label "wf_human_str"
    cpus 1
    memory 4.GB
    input:
        tuple val(join_key), path(vcf), path(tsv)
        path(variant_catalogue_hg38)
    output:
        tuple val(join_key), path("${join_key}_repeat-expansion_annotated.vcf.gz"), path("${join_key}_repeat-expansion_annotated.vcf.gz.tbi"), path("*_plot.tsv"), path("*_annotated.tsv")
    script:
        def chr = join_key
        """
        stranger -f ${variant_catalogue_hg38} ${vcf} \
            | sed 's/\\ /_/g' \
            | bgzip -c > ${chr}_repeat-expansion_annotated.vcf.gz
        tabix -p vcf ${chr}_repeat-expansion_annotated.vcf.gz
        SnpSift extractFields ${chr}_repeat-expansion_annotated.vcf.gz \
            CHROM POS ALT FILTER REF RL RU REPID VARID STR_STATUS > ${chr}_repeat-expansion_annotated.tsv
        SnpSift extractFields ${chr}_repeat-expansion_annotated.vcf.gz \
            CHROM POS DisplayRU STR_NORMAL_MAX STR_PATHOLOGIC_MIN VARID Disease > ${chr}_repeat-expansion_plot.tsv
        """
}


process bam_region_filter {
    // subset contig BAM to include only STR regions
    cpus 1
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        path (repeat_bed)
    output:
        tuple path("*str_regions.bam"), path("*str_regions.bam.bai"), val(sub_meta), emit: region_bam, optional: true
    script:
        sub_meta = [id: xam_meta.id, sq: xam_meta.sq, alias: xam_meta.alias]
        """
        { grep ${xam_meta.sq} -Fw ${repeat_bed} || true; } > repeats_subset.bed

        if [[ -s repeats_subset.bed ]]; then
            samtools view -b -h --write-index -o "${xam_meta.sq}.wf_str_regions.bam##idx##${xam_meta.sq}.wf_str_regions.bam.bai" -L repeats_subset.bed ${xam}
        else
            echo "blank subset BED"
        fi
        """
}

process bam_read_filter {
    // subset contig STR regions BAM to include only supporting reads from straglr
    cpus 1
    memory 4.GB
    input:
        tuple val(chr), path(xam), path(xam_idx), val(xam_meta), path(vcf), path(straglr_tsv)
    output:
        tuple path ("*str_reads.bam"), path("*str_reads.bam.bai"), val(xam_meta)
    shell:
        """
        tail -n +3 !{straglr_tsv} | cut -f6 > reads_to_filter.txt
        samtools view --write-index -N reads_to_filter.txt -o !{chr}.wf_str_reads.bam##idx##!{chr}.wf_str_reads.bam.bai !{xam}
        """
}

process generate_str_content {
    // extract content info from BAM and generate TSV files for plot data
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        tuple val(chr), path(straglr_vcf), path(straglr_tsv), path(annotated_vcf), path(annotated_vcf_tbi), path(stranger_tsv), path(stranger_annotation), path(xam), path(xam_idx), val(xam_meta)
        path (repeat_bed)
    output:
        tuple val(sub_meta), path ("*str-content.csv"), optional: true
    script:
        sub_meta = ["alias": xam_meta.alias]
        """
        workflow-glue generate_str_content \
            --straglr ${straglr_tsv} \
            --stranger ${stranger_tsv} \
            --chr ${chr} \
            --repeat_bed ${repeat_bed} \
            --str_reads_bam ${xam}
        """
}

process merge_tsv {
    // merge the contig TSVs/CSVs
    cpus 1
    memory 4.GB
    input:
        path (plot_tsv)
        path (straglr_tsv)
        path (stranger_tsv)
        tuple val(xam_meta), path(str_content_csv)
    output:
        tuple path ("*plot.tsv"), path ("*straglr.tsv"), path ("*stranger.tsv"), path ("*str-content-all.csv")
    script:
        """
        awk 'NR == 1 || FNR > 1' ${plot_tsv} >${xam_meta.alias}_plot.tsv
        awk 'NR == 1 || FNR > 1' ${stranger_tsv} >${xam_meta.alias}_stranger.tsv
        # Ignore the first line from straglr_tsv as it has a two-line header. Avoid sed altogether
        awk 'NR == 2 || FNR > 2' ${straglr_tsv} > ${xam_meta.alias}.wf_str.straglr.tsv
        awk 'NR == 1 || FNR > 1' ${str_content_csv} > ${xam_meta.alias}_str-content-all.csv
        """
}


process make_report {
    label "wf_common"
    cpus 1
    memory 16.GB
    input:
        tuple val(xam_meta), path(vcf), path(vcf_idx)
        path(straglr_tsv)
        path(plot_tsv)
        path(stranger_annotation)
        path(str_content)
        path "versions.txt"
        path "params.json"
        path(bam_stats)
        val(sex)
    output:
        path "*wf-human-str-report.html", emit: html
    script:
        def report_name = "${xam_meta.alias}.wf-human-str-report.html"
        // if params.sex is not provided, assume the workflow inferred it
        String sex_source = params.sex ? "user-provided" : "workflow-inferred"
        """
        workflow-glue report_str \
            -o $report_name \
            --params params.json \
            --sample_name ${xam_meta.alias} \
            --version versions.txt \
            --vcf ${vcf} \
            --straglr ${straglr_tsv} \
            --stranger ${plot_tsv} \
            --stranger_annotation ${stranger_annotation} \
            --str_content ${str_content} \
            --read_stats ${bam_stats} \
            --sex ${sex} \
            --sex_source ${sex_source} \
            --workflow_version ${workflow.manifest.version}
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
