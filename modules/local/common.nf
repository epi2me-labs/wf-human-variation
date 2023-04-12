import groovy.json.JsonBuilder

process cram_cache {
    input:
        path reference
    output:
        path("ref_cache/"), emit: cram_cache
    script:
    def is_conda = workflow.profile.toLowerCase().contains("conda")
    """
    if [[ "${is_conda}" == "true" ]]; then
        wget https://raw.githubusercontent.com/samtools/samtools/master/misc/seq_cache_populate.pl;
        # Invoke without messing with PATH and +x
        INVOCATION='perl seq_cache_populate.pl'
    else
        # Invoke from binary installed to container PATH
        INVOCATION='seq_cache_populate.pl'
    fi
    \$INVOCATION -root ref_cache/ ${reference}
    """
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
        tuple path(ref), path(ref_idx), path(ref_cache)
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
        # extract first 3 columns of input BED to prevent col 4 leaking as label into outputs [CW-1702]
        awk -v OFS='\\t' '{print \$1, \$2, \$3}' ${target_bed} > cut.bed
        export REF_PATH=${ref}
        export MOSDEPTH_PRECISION=3
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
        tuple path(ref), path(ref_idx), path(ref_cache)
    output:
        tuple \
            path("${params.sample_name}.mapula.csv"),
            path("${params.sample_name}.mapula.json")
    script:
        def ref_path = "${ref_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
        """
        export REF_PATH=${ref_path}
        mapula count -r ${ref} -f all -n '${params.sample_name}.mapula' ${xam}
        """
}


process readStats {
    label "wf_human_snp"
    cpus 4
    input:
        tuple path(xam), path(xam_idx), val(xam_meta)
        path target_bed
        tuple path(ref), path(ref_idx), path(ref_cache)
    output:
        path "${params.sample_name}.readstats.tsv.gz"
    script:
        def ref_path = "${ref_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
        """
        export REF_PATH="${ref_path}"
        bamstats --threads 3 "${xam}" | gzip > "${params.sample_name}.readstats.tsv.gz"
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
        tuple path(reference), path(ref_idx), path(ref_cache)
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    faidx --transform bed $reference > allChromosomes.bed
    """
}


//TODO reference is later copied to out_dir to save us a lot of trouble but is wasteful
//     additionally, --reference is hacked to allow the actual_ref to be opened
//     (as it cannot be guaranteed to exist in the out_dir at this point)
//TODO alignment track only output when alignment has been done, for now
//TODO --variant locations should be constructed legitimately instead of guessed
process configure_jbrowse {
    input:
        tuple path(reference), path(ref_idx), path(ref_cache)
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
        --reference ${reference} ${params.out_dir}/${reference.name} ${params.out_dir}/${ref_idx.name} \
        ${alignment} \
        ${snp_variant} \
        ${sv_variant} > jbrowse.json
    """
}


process minimap2_ubam {
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        path reference
        path old_reference
        tuple path(reads), path(reads_idx)
    output:
        tuple path("${params.sample_name}.cram"), path("${params.sample_name}.cram.crai"), emit: alignment
    script:
    def bam2fq_ref = old_reference.name != "OPTIONAL_FILE" ? "--reference ${old_reference}" : ''
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${bam2fq_ref} ${reads} | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${reference} - \
    | samtools sort -@ ${params.ubam_sort_threads} --write-index -o ${params.sample_name}.cram##idx##${params.sample_name}.cram.crai -O CRAM --reference ${reference} -
    """
}


process get_genome {
    cpus 1
    input:
        tuple path(bam), path(bai)
    output:
        path("${bam}_genome.txt")
        path("output.txt"), emit: genome_label optional true
        env genome_string, emit: genome_build optional true
     script:
        """
        samtools idxstats ${bam} > ${bam}_genome.txt
        get_genome.py --chr_counts ${bam}_genome.txt -o output.txt -w str
        genome_string=`cat output.txt`
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


// Make bam QC reporting.
process failedQCReport  {
    cpus 1
    input:
        tuple val(cvg), path(xam), path(xam_idx), val(xam_meta)
        file read_summary
        tuple path(mosdepth_bed), path(mosdepth_dist), path(mosdepth_threshold) // MOSDEPTH_TUPLE
        path versions
        path "params.json"
    output:
        path "*report.html"
    script:
        report_name = "${params.sample_name}.wf-human-qc-report.html"
        wfversion = workflow.manifest.version
        if( workflow.commitId ){
            wfversion = workflow.commitId
        }
        """
        workflow-glue report_qc \
            $report_name \
            --versions $versions \
            --params params.json \
            --read_stats $read_summary \
            --read_depth $mosdepth_dist \
            --revision $wfversion \
            --commit $workflow.commitId \
            --low_cov ${params.bam_min_coverage}
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
