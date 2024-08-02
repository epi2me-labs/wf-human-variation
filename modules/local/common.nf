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
        val (create_gene_summary)
    output:
        tuple \
            val(xam_meta),
            path("${xam_meta.alias}.regions.bed.gz"),
            path("${xam_meta.alias}.mosdepth.global.dist.txt"),
            path("${xam_meta.alias}.thresholds.bed.gz"), emit: mosdepth_tuple
        path "${xam_meta.alias}.mosdepth.summary.txt", emit: summary
        path("${xam_meta.alias}.per-base.bedgraph.gz"), emit: perbase, optional: true
        path "${xam_meta.alias}.gene_summary.tsv", emit: gene_summary, optional: true
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
        ${xam_meta.alias} \
        $xam

        # Rename the output, avoiding ambiguity in the output formatting
        if [ -e ${xam_meta.alias}.per-base.bed.gz ]; then
            mv ${xam_meta.alias}.per-base.bed.gz ${xam_meta.alias}.per-base.bedgraph.gz
        fi

        # If gene summary requested and a BED provided, run mosdepth again without -x to get precise coverage 
        # Use thresholds and regions file to create gene summary
        if [ "${create_gene_summary}" = true ]; then
            mosdepth \
            -t $task.cpus \
            -b ${target_bed} \
            --thresholds 1,10,15,20,30 \
            --no-per-base \
            ${xam_meta.alias}.gene \
            $xam

            gunzip -c ${xam_meta.alias}.gene.thresholds.bed.gz > thresholds.bed
            gunzip -c ${xam_meta.alias}.gene.regions.bed.gz  > regions.bed

            workflow-glue generate_gene_summary --mosdepth_threshold thresholds.bed --mosdepth_average regions.bed --output ${xam_meta.alias}.gene_summary.tsv
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
        path "${xam_meta.alias}.readstats.tsv.gz", emit: read_stats
        path "${xam_meta.alias}.flagstat.tsv", emit: flagstat
        path "${xam_meta.alias}-histograms/", emit: histograms
        path "${xam_meta.alias}.runids.txt", emit: runids
        path "${xam_meta.alias}.basecallers.txt", emit: basecallers
    script:
        def stats_threads = Math.max(task.cpus - 1, 1)
        """
        bamstats \
            -s ${xam_meta.alias} \
            -i ${xam_meta.alias}.per-file-runids.txt \
            -l ${xam_meta.alias}.basecallers.tsv \
            --histogram ${xam_meta.alias}-histograms \
            -u \
            -f ${xam_meta.alias}.flagstat.tsv \
            --threads ${stats_threads} \
            "${xam}" | gzip > "${xam_meta.alias}.readstats.tsv.gz"
        # get unique run IDs
        awk -F '\\t' '
            NR==1 {for (i=1; i<=NF; i++) {ix[\$i] = i}}
            # only print run_id if present
            NR>1 && \$ix["run_id"] != "" {print \$ix["run_id"]}
        ' ${xam_meta.alias}.per-file-runids.txt | sort | uniq > ${xam_meta.alias}.runids.txt
        # get unique basecall models
        awk -F '\\t' '
            NR==1 {for (i=1; i<=NF; i++) {ix[\$i] = i}}
            # only print basecall model if present
            NR>1 && \$ix["basecaller"] != "" {print \$ix["basecaller"]}
        ' ${xam_meta.alias}.basecallers.tsv | sort | uniq > ${xam_meta.alias}.basecallers.txt
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
    label "wf_common"
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
        tuple val(xam_meta),
            path(regions),
            path(dists),
            path(thresholds)

    output:
        path "${bed.baseName}.filt.bed", emit: filt_bed
        tuple \
            val(xam_meta),
            path("${xam_meta.alias}.regions.filt.bed.gz"),
            path("${xam_meta.alias}.mosdepth.global.dist.txt"),
            path("${xam_meta.alias}.thresholds.bed.gz"), emit: mosdepth_tuple
    shell:
    '''
    # Get intervals with average coverage above minimum
    zcat !{regions} | awk '$NF>=!{params.bam_min_coverage}' | bgzip -c > !{xam_meta.alias}.regions.filt.bed.gz
    
    # Extract original regions with reasonable coverage. We first intersect, sort the kept intervals,
    # merge the adjacent and then sort again.
    bedtools intersect -a !{bed} -b !{xam_meta.alias}.regions.filt.bed.gz | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - | \
        sort -k1,1 -k2,2n > !{bed.baseName}.filt.bed
    '''
}


// Make bam QC reporting.
process failedQCReport  {
    label "wf_common"
    cpus 1
    memory 12.GB
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats.tsv.gz"),
            path("flagstats.tsv"),
            path("hists"),
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
        def report_name = "${xam_meta.alias}.wf-human-alignment-report.html"
        """
        workflow-glue report_al \\
            --name ${report_name} \\
            --sample_name ${xam_meta.alias} \\
            --stats_fn readstats.tsv.gz \\
            --summary_dir summary/ \\
            --hists_dir hists/ \\
            --flagstat_fn flagstats.tsv \\
            --depths_dir depths/ \\
            --versions versions.txt \\
            --window_size ${params.depth_window_size} \\
            --params params.json \\
            ${genome_wide_depth} \\
            --low_cov ${params.bam_min_coverage} \\
            --workflow_version ${workflow.manifest.version}
        """
}

// Alignment report
process makeAlignmentReport {
    label "wf_common"
    cpus 1
    memory 12.GB
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats.tsv.gz"),
            path("flagstats.tsv"),
            path("hists"),
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
        def report_name = "${xam_meta.alias}.wf-human-alignment-report.html"
        """
        workflow-glue report_al \\
            --name ${report_name} \\
            --sample_name "${xam_meta.alias}" \\
            --stats_fn readstats.tsv.gz \\
            --summary_dir summary/ \\
            --hists_dir hists/ \\
            --reference_fai ref.fasta.fai \\
            --flagstat_fn flagstats.tsv \\
            --depths_dir depths/ \\
            --versions versions.txt \\
            --window_size ${params.depth_window_size} \\
            --params params.json \\
            --workflow_version ${workflow.manifest.version}
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
        tuple val(xam_meta), path("input.vcf.gz"), path("input.vcf.gz.tbi"), val(contig)
        val(genome)
        val(output_label)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_${output_label}*.vcf.gz"), path("${xam_meta.alias}.wf_${output_label}*.vcf.gz.tbi"), emit: annot_vcf
    shell:
    '''
    if [ "!{contig}" == '*' ]; then
        # SV is quick to annotate, dont bother breaking it apart
        INPUT_FILENAME=input.vcf.gz
        OUTPUT_LABEL="!{output_label}"
    else
        # SNP is slow to annotate, we'll break it apart by contig
        # and merge it back later. filter the master VCF to current contig
        # Also, relabel chrM as chrMT so SnpEff can annotate
        bcftools view -r !{contig} input.vcf.gz | sed 's/^chrM\\t/chrMT\\t/' | bgzip > input.chr.vcf.gz

        INPUT_FILENAME=input.chr.vcf.gz
        OUTPUT_LABEL="!{output_label}.!{contig}"
    fi

    # deal with samples which aren't hg19 or hg38
    if [[ "!{genome}" != "hg38" ]] && [[ "!{genome}" != "hg19" ]]; then
        # return the original VCF and index as the outputs
        cp ${INPUT_FILENAME} !{xam_meta.alias}.wf_${OUTPUT_LABEL}.vcf.gz
        cp ${INPUT_FILENAME}.tbi !{xam_meta.alias}.wf_${OUTPUT_LABEL}.vcf.gz.tbi
    else
        # do some annotation
        if [[ "!{genome}" == "hg38" ]]; then
            snpeff_db="GRCh38.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh38.vcf.gz"

        elif [[ "!{genome}" == "hg19" ]]; then
            snpeff_db="GRCh37.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh37.vcf.gz"
        fi

        # Revert back chrMT to chrM
        snpEff -Xmx!{task.memory.giga - 1}g ann -noStats -noLog $snpeff_db ${INPUT_FILENAME} | sed 's/^chrMT\\t/chrM\\t/' > !{xam_meta.alias}.intermediate.snpeff_annotated.vcf

        # Add ClinVar annotations
        SnpSift annotate $clinvar_vcf !{xam_meta.alias}.intermediate.snpeff_annotated.vcf | bgzip > !{xam_meta.alias}.wf_${OUTPUT_LABEL}.vcf.gz
        tabix !{xam_meta.alias}.wf_${OUTPUT_LABEL}.vcf.gz

        # tidy up
        rm !{xam_meta.alias}.intermediate*
    fi
    '''
}

process sift_clinvar_vcf {
    label "snpeff_annotation"
    cpus 1
    memory 3.GB
    input:
        tuple val(xam_meta), path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(genome)
        val(output_label)
    output:
        path("${xam_meta.alias}.wf_${output_label}_clinvar.vcf"), emit: final_vcf_clinvar
    shell:
    '''
    # deal with samples which aren't hg19 or hg38
    if [[ "!{genome}" != "hg38" ]] && [[ "!{genome}" != "hg19" ]]; then
        # create an empty ClinVar VCF
        touch !{xam_meta.alias}.wf_!{output_label}_clinvar.vcf
    else
        bcftools view input.vcf.gz | SnpSift filter "( exists CLNSIG )" > !{xam_meta.alias}.wf_!{output_label}_clinvar.vcf
    fi
    '''
}

process concat_vcfs {
    cpus 2
    memory 3.GB
    input:
        tuple val(xam_meta), path (vcfs_artifacts, stageAs: "vcfs/*")
        val(output_label)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.${output_label}.vcf.gz"), path("${xam_meta.alias}.${output_label}.vcf.gz.tbi"), emit: final_vcf
    script:
        def concat_threads = Math.max(task.cpus - 1, 1)
        def alias_output_label = "${xam_meta.alias}.${output_label}"
        """
        bcftools concat --threads ${concat_threads} -O u vcfs/*.vcf.gz | bcftools sort -O z - > "${alias_output_label}.vcf.gz"
        tabix -p vcf ${alias_output_label}.vcf.gz
        """
}

process haploblocks {
    // Extract phased blocks from a VCF file
    cpus 1
    memory 8.GB
    input:
        tuple val(xam_meta), path(phased_vcf), path(phased_tbi)
        val output_label
    output:
        path "${xam_meta.alias}.wf_${output_label}.haploblocks.gtf", emit: phase_blocks
    script:
        """
        # Prepare correct input file
        whatshap stats --gtf=${xam_meta.alias}.wf_${output_label}.haploblocks.gtf ${phased_vcf}
        """
}

process bed_filter {
    // filter a BED/VCF/GFF using a BED file
    cpus 1
    memory 4.GB
    input:
        tuple val(xam_meta), path(input, stageAs: 'input.gz'), path(input_tbi, stageAs: 'input.gz.tbi')
        path(bed)
        val(subworkflow)
        val(file_type)
    output:
        tuple val(xam_meta), path("${xam_meta.alias}.wf_${subworkflow}.${file_type}.gz"), path("${xam_meta.alias}.wf_${subworkflow}.${file_type}.gz.tbi"), emit: filtered
    script:
        """
        bedtools intersect -u -header -a input.gz -b ${bed} | bgzip -c > ${xam_meta.alias}.wf_${subworkflow}.${file_type}.gz
        tabix ${xam_meta.alias}.wf_${subworkflow}.${file_type}.gz
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
// NOTE The keys in here are rather sad to look at but form part of downstream
//   processes for several 3rd party providers and MUST NOT be fiddled with.
// Currently, the workflow's sample_name parameter is used to fill in
//   meta.sample_sheet.alias, this should be populated more fully by an
//   input sample_sheet in a not so distant future (and sample_name deprecated)
process combine_metrics_json {
    label "wf_common"
    cpus 1
    memory 4.GB
    input:
        path jsons
        path "flagstat.tsv"
        path "hists"
        tuple \
            val(xam_meta),
            path("regions.bed.gz"),
            path("mosdepth.global.dist.txt"),
            path("thresholds.bed.gz")
        path "mosdepth.summary.txt"
        path haplocheck, stageAs: "haplocheck/*"
        val(sex)
    output:
        path "${xam_meta.alias}.stats.json", emit: json
    script:
        String input_jsons = jsons.name != 'OPTIONAL_FILE' ? "--jsons ${jsons}" : ""
        // only emit an inferred sex if sex is defined and params.sex is not
        String sex_arg = (!params.sex && sex) ? "--inferred_sex ${sex}" : ""
        // If haplocheck is optional file, then skip it
        String haplocheck_arg = params.haplocheck && haplocheck.baseName != "OPTIONAL_FILE" ? "--haplocheck ${haplocheck}" : ""
        """
        workflow-glue combine_jsons \
            --bamstats_flagstats flagstat.tsv \
            --bamstats_hists hists \
            --mosdepth_summary mosdepth.summary.txt \
            --mosdepth_thresholds thresholds.bed.gz \
            ${haplocheck_arg} \
            ${sex_arg} \
            ${input_jsons} \
            --metadata "sample_sheet.alias=${xam_meta.alias}" \
            --output ${xam_meta.alias}.stats.json
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

process haplocheck {
    // Estimate contamination level from the mitogenome
    // using haplocheck.
    cpus 2
    memory 8.GB
    errorStrategy 'ignore'
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        val mt_chr
    output:
        path "${xam_meta.alias}.haplocheck.tsv"
    script:
        """
        # Create the reference genome
        samtools faidx ${ref} ${mt_chr} > mt.fa && samtools faidx mt.fa

        # Extract mito-genome reads first
        samtools view -@${task.cpus - 1} -hb ${xam} ${mt_chr} > ${xam_meta.alias}
        samtools index -@${task.cpus} ${xam_meta.alias}

        # Count how many reads are mapped to the mitogenome.
        n_reads=\$( samtools view -c ${xam_meta.alias} )

        # Run the commands if there are reads to process.
        # We already checked if the reference has a valid mitogenome.
        if [ \$n_reads -gt 0 ]; then
            # Run mutserve
            java -jar `which mutserve.jar` call \
                --level 0.01 \
                --reference mt.fa \
                --mapQ 20 \
                --baseQ 20 \
                --output mt.vcf.gz \
                --no-ansi \
                --threads ${task.cpus} \
                ${xam_meta.alias}

            # Before running haplocheck, count how many sites are in the
            # VCF file. If it is 0, then do not proceed.
            n_sites=0
            if [ -e mt.vcf.gz ]; then
                n_sites=\$( zcat mt.vcf.gz | awk '\$1!~"#" {print}' | wc -l | awk '{print \$1}' )
            fi

            # Run haplocheck if there are sites in the VCF file
            if [ \$n_sites -gt 0 ]; then
                haplocheck --out ${xam_meta.alias}.haplocheck.tsv mt.vcf.gz
            # Otherwise, save as NV (no value) as opposed to ND (not determined)
            else
                echo "Sample\tContamination Status\tContamination Level\tDistance\tSample Coverage" > ${xam_meta.alias}.haplocheck.tsv
                echo "${xam_meta.alias}\tNA\tNV\t0\t0" >> ${xam_meta.alias}.haplocheck.tsv
            fi
        # If no reads are found, create the placeholder.
        else
            echo "Sample\tContamination Status\tContamination Level\tDistance\tSample Coverage" > ${xam_meta.alias}.haplocheck.tsv
            echo "${xam_meta.alias}\tNA\tNV\t0\t0" >> ${xam_meta.alias}.haplocheck.tsv
        fi
        """
}
