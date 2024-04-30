import groovy.json.JsonBuilder

process cram_cache {
    cpus 1
    memory 4.GB
    input:
        path reference
    output:
        path("ref_cache/"), emit: ref_cache
        env(REF_PATH), emit: ref_path
    shell:
    '''
    # Invoke from binary installed to container PATH
    seq_cache_populate.pl -root ref_cache/ !{reference}
    REF_PATH="ref_cache/%2s/%2s/%s"
    '''
}

process index_ref_fai {
    cpus 1
    memory 4.GB
    input:
        file reference
    output:
        path "${reference}.fai", emit: reference_index
    """
    samtools faidx ${reference}
    """
}

process index_ref_gzi {
    cpus 1
    memory 4.GB
    input:
        file reference
    output:
        path "${reference}.gzi", emit: reference_index
    """
    bgzip -r ${reference}
    """
}

// NOTE -f required to compress symlink
process decompress_ref {
    cpus 1
    memory 4.GB
    input:
        file compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    """
    gzip -df ${compressed_ref}
    """
}


//NOTE grep MOSDEPTH_TUPLE if changing output tuple
process mosdepth {
    cpus 2
    memory {4.GB * task.attempt}
    maxRetries 2
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        file target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH) 
        val (window_size)
    output:
        tuple \
            path("${params.sample_name}.regions.bed.gz"),
            path("${params.sample_name}.mosdepth.global.dist.txt"),
            path("${params.sample_name}.thresholds.bed.gz"), emit: mosdepth_tuple
        path "${params.sample_name}.mosdepth.summary.txt", emit: summary
        path("${params.sample_name}.per-base.bedgraph.gz"), emit: perbase, optional: true
    script:
        def perbase_args = params.depth_intervals ? "" : "--no-per-base"
        """
        export REF_PATH=${ref}
        export MOSDEPTH_PRECISION=3
        # extract first 3 columns of input BED to prevent col 4 leaking as label into outputs [CW-1702]
        # and convert them into windows of the given size [CW-2015]
        # The workflow now sort the bed input, merge overlapping intervals and then build windows
        # preventing crash in downstream tools [CW-2247]
        sort -k 1,1 -k2,2n ${target_bed} | \
            bedtools merge -i - | \
            bedtools makewindows -b - -w ${window_size} > cut.bed
        # Run mosdepth
        mosdepth \
        -x \
        -t $task.cpus \
        -b cut.bed \
        --thresholds 1,10,15,20,30 \
        ${perbase_args} \
        ${params.sample_name} \
        $xam

        # Rename the output, avoiding ambiguity in the output formatting
        if [ -e ${params.sample_name}.per-base.bed.gz ]; then
            mv ${params.sample_name}.per-base.bed.gz ${params.sample_name}.per-base.bedgraph.gz
        fi
        """
}


process readStats {
    label "wf_common"
    cpus 4
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        path target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "${params.sample_name}.readstats.tsv.gz", emit: read_stats
        path "${params.sample_name}.flagstat.tsv", emit: flagstat
    script:
        def stats_threads = Math.max(task.cpus - 1, 1)
        """
        bamstats -s ${params.sample_name} -u -f ${params.sample_name}.flagstat.tsv --threads ${stats_threads} "${xam}" | gzip > "${params.sample_name}.readstats.tsv.gz"
        """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish_artifact {
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


process getAllChromosomesBed {
    cpus 1
    memory 4.GB
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    shell:
    '''
    awk '{OFS="\t"; print $1, "0", $2}' !{ref_idx} > allChromosomes.bed
    '''
}


//TODO reference is later copied to out_dir to save us a lot of trouble but is wasteful
//     additionally, --reference is hacked to allow the actual_ref to be opened
//     (as it cannot be guaranteed to exist in the out_dir at this point)
//TODO alignment track only output when alignment has been done, for now
//TODO --variant locations should be constructed legitimately instead of guessed
process configure_jbrowse {
    cpus 2
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        path("jbrowse.json")
    script:
    boolean output_bam = xam_meta.output
    def snp_variant = params.snp ? "--variant snp '${params.out_dir}/${params.sample_name}.wf_snp.vcf.gz' '${params.out_dir}/${params.sample_name}.wf_snp.vcf.gz.tbi'" : ''
    def sv_variant = params.sv ? "--variant sv '${params.out_dir}/${params.sample_name}.wf_sv.vcf.gz' '${params.out_dir}/${params.sample_name}.wf_sv.vcf.gz.tbi'" : ''
    def alignment = output_bam ? "--alignment '${params.out_dir}/${xam.name}' '${params.out_dir}/${xam_idx.name}'" : ''
    """
    workflow-glue configure_jbrowse \
        --reference ${ref} '${params.out_dir}/${ref.name}' '${params.out_dir}/${ref_idx.name}' \
        ${alignment} \
        ${snp_variant} \
        ${sv_variant} > jbrowse.json
    """
}

process getGenome {
    cpus 1
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        env genome_build, emit: genome_build optional true
     script:
        // set flags for subworkflows that have genome build restrictions
        def str_arg = params.str ? "--str" : ""
        def cnv_arg = params.cnv ? "--cnv" : ""
        def qdnaseq_arg = params.use_qdnaseq ? "--use_qdnaseq" : ""
        """
        # use view -H rather than idxstats, as idxstats will still cause a scan of the whole CRAM (https://github.com/samtools/samtools/issues/303)
        samtools view -H ${xam} --no-PG | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*LN:([^[:space:]]*).*,\\1\\t\\2,p' > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt ${str_arg} ${cnv_arg} ${qdnaseq_arg}
        genome_build=`cat output.txt`
        """
}

process eval_downsampling {
    cpus 1
    memory 4.GB
    input:
        path mosdepth_summary
        path bed
    output:
        path "ratio.txt", emit: downsampling_ratio
    script:
        def with_bed = bed.name != 'OPTIONAL_FILE' ? "--bed ${bed}" : ""
        """
        workflow-glue downsampling_ratio \
            --downsample_depth ${params.downsample_coverage_target} \
            --margin ${params.downsample_coverage_margin} \
            --summary ${mosdepth_summary} \
            ${with_bed} > ratio.txt
        """
}


process downsampling {
    cpus 4
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(to_downsample), val(downsampling_rate)
        tuple val(xam_fmt), val(xai_fmt)
    output:
        tuple path("downsampled.${xam_fmt}"), path("downsampled.${xam_fmt}.${xai_fmt}"), val(xam_meta), emit: xam optional true
    script:
        """
        samtools view \\
            -@ ${task.cpus} \\
            -h \\
            -O ${xam_fmt} \\
            --write-index \\
            --reference $ref \\
            --subsample ${downsampling_rate} \\
            -o downsampled.${xam_fmt}##idx##downsampled.${xam_fmt}.${xai_fmt} \\
            $xam
        """
}


// Process to get the genome coverage from the mosdepth summary.
process get_region_coverage {
    cpus 1
    memory 4.GB
    input:
        path bed
        tuple path(regions),
            path(dists),
            path(thresholds)

    output:
        path "${bed.baseName}.filt.bed", emit: filt_bed
        tuple \
            path("${params.sample_name}.regions.filt.bed.gz"),
            path("${params.sample_name}.mosdepth.global.dist.txt"),
            path("${params.sample_name}.thresholds.bed.gz"), emit: mosdepth_tuple
    shell:
    '''
    # Get intervals with average coverage above minimum
    zcat !{regions} | awk '$NF>=!{params.bam_min_coverage}' | bgzip -c > !{params.sample_name}.regions.filt.bed.gz
    
    # Extract original regions with reasonable coverage. We first intersect, sort the kept intervals,
    # merge the adjacent and then sort again.
    bedtools intersect -a !{bed} -b !{params.sample_name}.regions.filt.bed.gz | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - | \
        sort -k1,1 -k2,2n > !{bed.baseName}.filt.bed
    '''
}


// Make bam QC reporting.
process failedQCReport  {
    cpus 1
    memory 12.GB
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats.tsv.gz"),
            path("flagstats/*"),
            path("depths/*"),
            path("summary/*"),
            path('ref.fasta'),
            path('ref.fasta.fai'),
            path('ref_cache/'),
            env(REF_PATH),
            path("versions.txt"),
            path("params.json")

    output:
        path "*.html"

    script:
        // CW-2585: When we provide a fai file, the workflow will add the missing intervals to reach the end of the chromosome.
        // If we do not provide that, then the interval will be considered as 0>BP, where BP is the highest value in the bed 
        // file for each chromosome. This will display the intervals not in the context of the chromosome (so showing a peak in
        // a small region, and flat everywhere else) but only for the regions selected.
        def genome_wide_depth = params.bed ? "" : "--reference_fai ref.fasta.fai"
        def report_name = "${params.sample_name}.wf-human-alignment-report.html"
        """
        workflow-glue report_al \\
            --name ${report_name} \\
            --sample_name ${params.sample_name} \\
            --stats_fn readstats.tsv.gz \\
            --summary_dir summary/ \\
            --flagstat_dir flagstats/ \\
            --depths_dir depths/ \\
            --versions versions.txt \\
            --window_size ${params.depth_window_size} \\
            --params params.json \\
            ${genome_wide_depth} \\
            --low_cov ${params.bam_min_coverage}
        """
}

// Alignment report
process makeAlignmentReport {
    cpus 1
    memory 12.GB
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats.tsv.gz"),
            path("flagstats/*"),
            path("depths/*"),
            path("summary/*"),
            path('ref.fasta'),
            path('ref.fasta.fai'),
            path('ref_cache/'),
            env(REF_PATH),
            path("versions.txt"),
            path("params.json")

    output:
        path "*.html"

    script:
        def report_name = "${params.sample_name}.wf-human-alignment-report.html"
        """
        workflow-glue report_al \\
            --name ${report_name} \\
            --sample_name "${params.sample_name}" \\
            --stats_fn readstats.tsv.gz \\
            --summary_dir summary/ \\
            --reference_fai ref.fasta.fai \\
            --flagstat_dir flagstats/ \\
            --depths_dir depths/ \\
            --versions versions.txt \\
            --window_size ${params.depth_window_size} \\
            --params params.json 
        """
}



// Create version of bamstats and mosdepth, as
// well as get parameters for reporting.
process getVersions {
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        mosdepth --version | sed 's/ /,/' >> versions.txt
        bamstats --version | awk '{print "bamstats,"\$1}' >> versions.txt
        """
}


process getParams {
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


process annotate_vcf {
    // use SnpEff to generate basic functional annotations,
    // followed by SnpSift annotate to add ClinVar annotations
    label "snpeff_annotation"
    cpus 1
    memory 6.GB
    input:
        tuple path("input.vcf.gz"), path("input.vcf.gz.tbi"), val(contig)
        val(genome)
        val(output_label)
    output:
        tuple path("${params.sample_name}.wf_${output_label}*.vcf.gz"), path("${params.sample_name}.wf_${output_label}*.vcf.gz.tbi"), emit: annot_vcf
    shell:
    '''
    if [ "!{contig}" == '*' ]; then
        # SV is quick to annotate, dont bother breaking it apart
        INPUT_FILENAME=input.vcf.gz
        OUTPUT_LABEL="!{output_label}"
    else
        # SNP is slow to annotate, we'll break it apart by contig
        # and merge it back later. filter the master VCF to current contig
        bcftools view -r !{contig} input.vcf.gz | bgzip > input.chr.vcf.gz
        INPUT_FILENAME=input.chr.vcf.gz
        OUTPUT_LABEL="!{output_label}.!{contig}"
    fi

    # deal with samples which aren't hg19 or hg38
    if [[ "!{genome}" != "hg38" ]] && [[ "!{genome}" != "hg19" ]]; then
        # return the original VCF and index as the outputs
        cp ${INPUT_FILENAME} !{params.sample_name}.wf_${OUTPUT_LABEL}.vcf.gz
        cp ${INPUT_FILENAME}.tbi !{params.sample_name}.wf_${OUTPUT_LABEL}.vcf.gz.tbi
    else
        # do some annotation
        if [[ "!{genome}" == "hg38" ]]; then
            snpeff_db="GRCh38.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh38.vcf.gz"

        elif [[ "!{genome}" == "hg19" ]]; then
            snpeff_db="GRCh37.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh37.vcf.gz"
        fi

        snpEff -Xmx!{task.memory.giga - 1}g ann -noStats -noLog $snpeff_db ${INPUT_FILENAME} > !{params.sample_name}.intermediate.snpeff_annotated.vcf
        # Add ClinVar annotations
        SnpSift annotate $clinvar_vcf !{params.sample_name}.intermediate.snpeff_annotated.vcf | bgzip > !{params.sample_name}.wf_${OUTPUT_LABEL}.vcf.gz
        tabix !{params.sample_name}.wf_${OUTPUT_LABEL}.vcf.gz

        # tidy up
        rm !{params.sample_name}.intermediate*
    fi
    '''
}

process sift_clinvar_vcf {
    label "snpeff_annotation"
    cpus 1
    memory 3.GB
    input:
        tuple path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(genome)
        val(output_label)
    output:
        path("${params.sample_name}.wf_${output_label}_clinvar.vcf"), emit: final_vcf_clinvar
    shell:
    '''
    # deal with samples which aren't hg19 or hg38
    if [[ "!{genome}" != "hg38" ]] && [[ "!{genome}" != "hg19" ]]; then
        # create an empty ClinVar VCF
        touch !{params.sample_name}.wf_!{output_label}_clinvar.vcf
    else
        bcftools view input.vcf.gz | SnpSift filter "( exists CLNSIG )" > !{params.sample_name}.wf_!{output_label}_clinvar.vcf
    fi
    '''
}

process concat_vcfs {
    cpus 2
    memory 3.GB
    input:
        path (vcfs_artifacts, stageAs: "vcfs/*")
        val(prefix)
    output:
        tuple path ("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: final_vcf
    script:
        def concat_threads = Math.max(task.cpus - 1, 1)
        """
        bcftools concat --threads ${concat_threads} -O u vcfs/*.vcf.gz | bcftools sort -O z - > "${prefix}.vcf.gz"
        tabix -p vcf ${prefix}.vcf.gz
        """
}

process haploblocks {
    // Extract phased blocks from a VCF file
    cpus 1
    memory 8.GB
    input:
        tuple path(phased_vcf), path(phased_tbi)
        val output_label
    output:
        path "${params.sample_name}.wf_${output_label}.haploblocks.gtf", emit: phase_blocks
    script:
        """
        # Prepare correct input file
        whatshap stats --gtf=${params.sample_name}.wf_${output_label}.haploblocks.gtf ${phased_vcf}
        """
}

process bed_filter {
    // filter a BED/VCF/GFF using a BED file
    cpus 1
    memory 4.GB
    input:
        tuple path(input, stageAs: 'input.gz'), path(input_tbi, stageAs: 'input.gz.tbi')
        path(bed)
        val(subworkflow)
        val(file_type)
    output:
        tuple path("${params.sample_name}.wf_${subworkflow}.${file_type}.gz"), path("${params.sample_name}.wf_${subworkflow}.${file_type}.gz.tbi"), emit: filtered
    script:
        """
        bedtools intersect -u -header -a input.gz -b ${bed} | bgzip -c > ${params.sample_name}.wf_${subworkflow}.${file_type}.gz
        tabix ${params.sample_name}.wf_${subworkflow}.${file_type}.gz
        """
}

process sanitise_bed {
    // Sanitise a BED file:
    // 1) check that chromosome labelling matches reference, exit if it doesn't
    // 1) check and fix whitespace
    // 2) sort by chromosome/position
    cpus 1
    memory 4.GB
    input:
        path(bed)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "${bed.baseName}.sanitised.bed"
    script:
        """
        sanitise_bed.py --ref ${ref} --bed ${bed} --bed_out ${bed.baseName}.sanitised.bed
        """
}

// Combine the JSON with base metrics for read stats,
// coverage, SNPs and SVs
process combine_metrics_json {
    cpus 1
    memory 4.GB
    input:
        path jsons
        path "readstats.tsv.gz"
        path "flagstat.tsv"
        tuple \
            path("regions.bed.gz"),
            path("mosdepth.global.dist.txt"),
            path("thresholds.bed.gz")
        path "mosdepth.summary.txt"
        val(sex)
    output:
        path "${params.sample_name}.stats.json", emit: json
    script:
        String input_jsons = jsons.name != 'OPTIONAL_FILE' ? "--jsons ${jsons}" : ""
        // only emit an inferred sex if sex is defined and params.sex is not
        String sex_arg = (!params.sex && sex) ? "--inferred_sex ${sex}" : ""
        """
        workflow-glue combine_jsons \
            --bamstats_readstats readstats.tsv.gz \
            --bamstats_flagstats flagstat.tsv \
            --mosdepth_summary mosdepth.summary.txt \
            --mosdepth_thresholds thresholds.bed.gz \
            ${sex_arg} \
            ${input_jsons} \
            --output ${params.sample_name}.stats.json
        """
}

// CNV output
// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_cnv {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}

process infer_sex {
    // Use relative coverage of chrX and chrY to infer sex
    //   inferred sex will be 'XY', 'XX'.
    // the workflow will err on calling XX if inference cannot be made
    //   to match the expected behaviour of the workflow before v2.1
    // The process will check if the coverage of the X chromosome
    // is substantially higher than the Y chromosome. We set
    // this value to 4 (i.e. cov of X 4x higer than cov of Y),
    // to account for the presence of pseudo-X region being covered
    // in the Y-chromosome.
    cpus 1
    input:
        path "mosdepth.summary.txt"
    output:
        env inferred_sex
    script:
        """
        # First, grep the X and Y coverage for the regions; then compute the rate if both are non-0; finally, infer the sex.
        inferred_sex=\$( 
            awk '
                BEGIN{x_cov=0; y_cov=0};
                \$1=="X_region" || \$1=="chrX_region" {x_cov=\$4};
                \$1=="Y_region" || \$1=="chrY_region" {y_cov=\$4};
                END {
                    if (x_cov > 0 && y_cov > 0) print x_cov/y_cov; else print 0
                }' mosdepth.summary.txt \
                | awk '\$1 > 4 {print "XX"}; \$1<=4 && \$1!=0 {print "XY"}; \$1==0 {print "XX"}' )
        """
}
