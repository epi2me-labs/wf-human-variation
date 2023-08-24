import groovy.json.JsonBuilder

process cram_cache {
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
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        file target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH) 
    output:
        tuple \
            path("${params.sample_name}.regions.bed.gz"),
            path("${params.sample_name}.mosdepth.global.dist.txt"),
            path("${params.sample_name}.thresholds.bed.gz"), emit: mosdepth_tuple
        path "${params.sample_name}.mosdepth.summary.txt", emit: summary
        path("${params.sample_name}.per-base.bed.gz"), emit: perbase, optional: true
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
            bedtools makewindows -b - -w ${params.depth_window_size} > cut.bed
        # Run mosdepth
        mosdepth \
        -x \
        -t $task.cpus \
        -b cut.bed \
        --thresholds 1,10,20,30 \
        ${perbase_args} \
        ${params.sample_name} \
        $xam
        """
}


process mapula {
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        path target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple \
            path("${params.sample_name}.mapula.csv"),
            path("${params.sample_name}.mapula.json")
    script:
        """
        mapula count -r ${ref} -f all -n '${params.sample_name}.mapula' ${xam}
        """
}


process readStats {
    label "wf_common"
    cpus 4
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        path target_bed
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "${params.sample_name}.readstats.tsv.gz", emit: read_stats
        path "${params.sample_name}.flagstat.tsv", emit: flagstat
    script:
        """
        bamstats -s ${params.sample_name} -u -f ${params.sample_name}.flagstat.tsv --threads 3 "${xam}" | gzip > "${params.sample_name}.readstats.tsv.gz"
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


// todo https://github.com/mdshw5/pyfaidx/pull/164
process getAllChromosomesBed {
    cpus 1
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    faidx --transform bed $ref > allChromosomes.bed
    """
}


//TODO reference is later copied to out_dir to save us a lot of trouble but is wasteful
//     additionally, --reference is hacked to allow the actual_ref to be opened
//     (as it cannot be guaranteed to exist in the out_dir at this point)
//TODO alignment track only output when alignment has been done, for now
//TODO --variant locations should be constructed legitimately instead of guessed
process configure_jbrowse {
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        path("jbrowse.json")
    script:
    boolean output_bam = xam_meta.output
    def snp_variant = params.snp ? "--variant snp ${params.out_dir}/${params.sample_name}.wf_snp.vcf.gz ${params.out_dir}/${params.sample_name}.wf_snp.vcf.gz.tbi" : ''
    def sv_variant = params.sv ? "--variant sv ${params.out_dir}/${params.sample_name}.wf_sv.vcf.gz ${params.out_dir}/${params.sample_name}.wf_sv.vcf.gz.tbi" : ''
    def alignment = output_bam ? "--alignment ${params.out_dir}/${xam.name} ${params.out_dir}/${xam_idx.name}" : ''
    """
    workflow-glue configure_jbrowse \
        --reference ${ref} "${params.out_dir}/${ref.name}" "${params.out_dir}/${ref_idx.name}" \
        ${alignment} \
        ${snp_variant} \
        ${sv_variant} > jbrowse.json
    """
}

process getGenome {
    cpus 1
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
    output:
        env genome_build, emit: genome_build optional true
     script:
        def workflow_arg = params.str ? "-w str" : ""
        """
        samtools idxstats ${xam} > ${xam}_genome.txt
        get_genome.py --chr_counts ${xam}_genome.txt -o output.txt ${workflow_arg}
        genome_build=`cat output.txt`
        """
}


process eval_downsampling {
    cpus 1
    input:
        path mosdepth_summary
    output:
        path "ratio.txt", emit: downsampling_ratio
    shell:
        '''
        workflow-glue downsampling_ratio \
            --downsample_depth !{params.downsample_coverage_target} \
            --margin !{params.downsample_coverage_margin} \
            --summary !{mosdepth_summary} > ratio.txt
        '''
}


process downsampling {
    cpus 4
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(to_downsample), val(downsampling_rate)
    output:
        tuple path("downsampled.bam"), path("downsampled.bam.bai"), val(xam_meta), emit: bam optional true
    script:
        """
        samtools view \\
            -@ ${task.cpus} \\
            -hb \\
            --subsample ${downsampling_rate} \\
            $xam > downsampled.bam
        samtools index -@ ${task.cpus} downsampled.bam
        """
}


// Process to get the genome coverage from the mosdepth summary.
process get_coverage {
    cpus 1
    input:
        path mosdepth_summary

    output:
        tuple env(passes), env(value), emit: pass
    shell:
        '''
        passes=$( awk 'BEGIN{v="false"}; NR>1 && $1=="total_region" && $4>=!{params.bam_min_coverage} {v="true"}; END {print v}' !{mosdepth_summary} )
        value=$( awk 'BEGIN{v=0}; NR>1 && $1=="total_region" {v=$4}; END {print v}' !{mosdepth_summary} )
        '''
}


// Process to get the genome coverage from the mosdepth summary.
process get_region_coverage {
    cpus 1
    input:
        path bed
        tuple path(regions),
            path(dists),
            path(thresholds)

    output:
        tuple env(passes), env(value), emit: pass
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

    # Return true if there are filtered intervals, otherwise false
    passes=$( zcat !{params.sample_name}.regions.filt.bed.gz | wc -l | awk 'BEGIN{v="false"}; $1>0 {v="true"}; END {print v}' )
    # If there are intervals, return average coverage, otherwise return 0
    value=$( zcat !{params.sample_name}.regions.filt.bed.gz | awk 'BEGIN{v=0; n=0}; {v+=$4; n+=1}; END {print v, n}' | awk '$2>0 {print $1/$2}; $2==0 {print 0}' )
    '''
}


// Make bam QC reporting.
process failedQCReport  {
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats/*"),
            path("flagstats/*"),
            path("depths/*"),
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
        """
        workflow-glue report_al \\
            --name wf-human-variation \\
            --stats_dir readstats/ \\
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
    input: 
        tuple path(xam),
            path(xam_idx),
            val(xam_meta),
            path("readstats/*"),
            path("flagstats/*"),
            path("depths/*"),
            path('ref.fasta'),
            path('ref.fasta.fai'),
            path('ref_cache/'),
            env(REF_PATH),
            path("versions.txt"),
            path("params.json")

    output:
        path "*.html"

    script:
        """
        workflow-glue report_al \\
            --name wf-human-variation \\
            --stats_dir readstats/ \\
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
    // use SnpEff to generate basic functional annotations, SnpSift annotate to add 
    // ClinVar annotations, and SnpSift filter to produce a separate VCF of Clinvar-annotated 
    // variants - if any variants are present in this file, it is used to populate a table in 
    // the report.
    label "snpeff_annotation"
    cpus {params.annotation_threads}
    input:
        tuple path("input.vcf.gz"), path("input.vcf.gz.tbi")
        val(genome)
        val(output_label)
    output:
        tuple path("${params.sample_name}.wf_${output_label}.vcf.gz"), path("${params.sample_name}.wf_${output_label}.vcf.gz.tbi"), emit: final_vcf
        path("${params.sample_name}.wf_${output_label}.snpEff_genes.txt")
        path("${params.sample_name}.wf_${output_label}_clinvar.vcf"), emit: final_vcf_clinvar
    shell:
    '''
    # deal with samples which aren't hg19 or hg38
    if [[ "!{genome}" != "hg38" ]] && [[ "!{genome}" != "hg19" ]]; then
        # return the original VCF and index as the outputs
        cp input.vcf.gz !{params.sample_name}.wf_!{output_label}.vcf.gz
        cp input.vcf.gz.tbi !{params.sample_name}.wf_!{output_label}.vcf.gz.tbi
        # create an empty snpEff_genes file
        touch !{params.sample_name}.wf_!{output_label}.snpEff_genes.txt
        # create an empty ClinVar VCF
        touch !{params.sample_name}.wf_!{output_label}_clinvar.vcf
    else
        # do some annotation
        if [[ "!{genome}" == "hg38" ]]; then
            snpeff_db="GRCh38.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh38.vcf.gz"
            
        elif [[ "!{genome}" == "hg19" ]]; then
            snpeff_db="GRCh37.p13"
            clinvar_vcf="${CLINVAR_PATH}/clinvar_GRCh37.vcf.gz"
        fi

        # Specify 4G of memory otherwise SnpEff will crash with the default 1G
        snpEff -Xmx4g ann $snpeff_db input.vcf.gz > !{params.sample_name}.snpeff_annotated.vcf
        # Add ClinVar annotations
        SnpSift annotate $clinvar_vcf !{params.sample_name}.snpeff_annotated.vcf > !{params.sample_name}.wf_!{output_label}.vcf
        # Get the ClinVar-annotated variants into a separate VCF
        cat !{params.sample_name}.wf_!{output_label}.vcf | SnpSift filter "( exists CLNSIG )" > !{params.sample_name}.wf_!{output_label}_clinvar.vcf
    
        bgzip -c !{params.sample_name}.wf_!{output_label}.vcf > !{params.sample_name}.wf_!{output_label}.vcf.gz
        tabix !{params.sample_name}.wf_!{output_label}.vcf.gz
    
        # tidy up
        rm !{params.sample_name}.snpeff_annotated.vcf
        rm !{params.sample_name}.wf_!{output_label}.vcf
        mv snpEff_genes.txt !{params.sample_name}.wf_!{output_label}.snpEff_genes.txt
    fi
    '''
}

process haploblocks {
    // Extract phased blocks from a VCF file
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
