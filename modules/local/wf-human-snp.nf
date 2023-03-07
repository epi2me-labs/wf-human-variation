import groovy.json.JsonBuilder


process make_chunks {
    // Do some preliminaries. Ordinarily this would setup a working directory
    // that all other commands would make use off, but all we need here are the
    // list of contigs and chunks.
    label "wf_human_snp"
    cpus 1
    input:
        tuple path(bam), path(bai)
        tuple path(ref), path(fai), path(ref_cache)
        path bed
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
    script:
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed}" : ''
        """
        mkdir -p clair_output
        python \$(which clair3.py) CheckEnvs \
            --bam_fn ${bam} \
            ${bedargs} \
            --output_fn_prefix clair_output \
            --ref_fn ${ref} \
            --vcf_fn ${params.vcf_fn} \
            --ctg_name ${params.ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs ${params.include_all_ctgs} \
            --threads 1  \
            --qual 2 \
            --sampleName ${params.sample_name} \
            --var_pct_full ${params.var_pct_full} \
            --ref_pct_full ${params.ref_pct_full} \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --min_contig_size ${params.min_contig_size}
        """
}


process pileup_variants {
    // Calls variants per region ("chunk") using pileup network.
    label "wf_human_snp"
    cpus 1
    errorStrategy 'retry'
    input:
        each region
        tuple path(bam), path(bai)
        tuple path(ref), path(fai), path(ref_cache)
        path model
    output:
        // TODO: make this explicit, why is pileup VCF optional?
        path "pileup_*.vcf", optional: true, emit: pileup_vcf_chunks
        path "gvcf_tmp_path/*", optional: true, emit: pileup_gvcf_chunks
    shell:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        // note: snp_min_af and indel_min_af have an impact on performance
        // TODO Fix REF_PATH
        '''
        export REF_PATH=!{ref_cache}/%2s/%2s/%s
        python $(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn !{model}/pileup \
            --bam_fn !{bam} \
            --call_fn pileup_!{region.contig}_!{region.chunk_id}.vcf \
            --ref_fn !{ref} \
            --ctgName !{region.contig} \
            --chunk_id !{region.chunk_id} \
            --chunk_num !{region.total_chunks} \
            --platform ont \
            --fast_mode False \
            --snp_min_af !{params.snp_min_af} \
            --indel_min_af !{params.indel_min_af} \
            --minMQ !{params.min_mq} \
            --minCoverage !{params.min_cov} \
            --call_snp_only False \
            --gvcf !{params.GVCF} \
            --temp_file_dir gvcf_tmp_path \
            --pileup
        ''' 
}


process aggregate_pileup_variants {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing.
    label "wf_human_snp"
    cpus 2
    input:
        tuple path(ref), path(fai), path(ref_cache)
        // these need to be named as original, as program uses info from
        // contigs file to filter
        path "input_vcfs/*"
        path contigs
    output:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi"), emit: pileup_vcf
        path "phase_qual", emit: phase_qual
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir input_vcfs/ \
            --vcf_fn_prefix pileup \
            --output_fn pileup.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        # Replaced bgzip with the faster bcftools index -n
        if [ "$( bcftools index -n pileup.vcf.gz )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

        bgzip -@ !{task.cpus} -fdc pileup.vcf.gz | \
            pypy $(which clair3.py) SelectQual --phase --output_fn .
        '''
}


process select_het_snps {
    // Filters a VCF by contig, selecting only het SNPs.
    label "wf_human_snp"
    cpus 2
    input:
        each contig
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectHetSnp.py#L29
        path "split_folder/phase_qual"
    output:
        tuple val(contig), path("split_folder/${contig}.vcf.gz"), path("split_folder/${contig}.vcf.gz.tbi"), emit: het_snps_vcf
    shell:
        '''
        pypy $(which clair3.py) SelectHetSnp \
            --vcf_fn pileup.vcf.gz \
            --split_folder split_folder \
            --ctgName !{contig}

        bgzip -c split_folder/!{contig}.vcf > split_folder/!{contig}.vcf.gz
        tabix split_folder/!{contig}.vcf.gz
        '''
}

process phase_contig_haplotag {
    // Tags reads in an input BAM from heterozygous SNPs
    // Also haplotag for those modes that need it (--str)

    label "wf_human_snp"
    cpus 4
    input:
        tuple val(contig), path(het_snps), path(het_snps_tbi), path(bam), path(bai), path(ref), path(fai), path(ref_cache)
    output:
        tuple val(contig), path("${contig}_hp.bam"), path("${contig}_hp.bam.bai"), path("phased_${contig}.vcf.gz"), emit: phased_bam_and_vcf
    shell:
        '''
            echo "Using whatshap for phasing"
            whatshap phase \
                --output phased_!{contig}.vcf.gz \
                --reference !{ref} \
                --chromosome !{contig} \
                --distrust-genotypes \
                --ignore-read-groups \
                !{het_snps} \
                !{bam}

        tabix -f -p vcf phased_!{contig}.vcf.gz

        whatshap haplotag \
            --reference !{ref} \
            --ignore-read-groups \
            --regions !{contig} \
            phased_!{contig}.vcf.gz \
            !{bam} \
        | samtools view -b -1 -@3 -o !{contig}_hp.bam

        # --write-index produces .csi not .bai, which downstream things seem not to like
        samtools index -@!{task.cpus} !{contig}_hp.bam

        '''
}

process phase_contig {
    // Tags reads in an input BAM from heterozygous SNPs
    // The haplotag step was removed in clair-v0.1.11 so this step re-emits
    //   the original BAM and BAI as phased_bam for compatability,
    //   but adds the VCF as it is now tagged with phasing information
    //   used later in the full-alignment model
    label "wf_human_snp"
    cpus 4
    input:
        tuple val(contig), path(het_snps), path(het_snps_tbi), path(bam), path(bai), path(ref), path(fai), path(ref_cache)
    output:
        tuple val(contig), path(bam), path(bai), path("phased_${contig}.vcf.gz"), emit: phased_bam_and_vcf
    shell:
        '''
        if [[ "!{params.use_longphase_intermediate}" == "true" ]]; then
            echo "Using longphase for phasing"
            # longphase needs decompressed 
            bgzip -@ !{task.cpus} -dc !{het_snps} > snps.vcf
            longphase phase --ont -o phased_!{contig} \
                -s snps.vcf -b !{bam} -r !{ref} -t !{task.cpus}
            bgzip phased_!{contig}.vcf
        else
            echo "Using whatshap for phasing"
            whatshap phase \
                --output phased_!{contig}.vcf.gz \
                --reference !{ref} \
                --chromosome !{contig} \
                --distrust-genotypes \
                --ignore-read-groups \
                !{het_snps} \
                !{bam}
        fi

        tabix -f -p vcf phased_!{contig}.vcf.gz
        '''
}

process merge_haplotagged_contigs {
    cpus params.threads
    // merge the haplotagged contigs to produce a single BAM file
    label "wf_human_snp"
    input:
        path contig_bams
    output:
        tuple path("*haplotagged.bam"), path("*haplotagged.bam.bai"), emit: merged_bam
    """
    samtools merge -@ $task.cpus -o ${params.sample_name}.haplotagged.bam ${contig_bams}
    samtools index -@ $task.cpus ${params.sample_name}.haplotagged.bam
    """
}

process get_qual_filter {
    // Determines quality filter for selecting candidate variants for second
    // stage "full alignment" calling.
    label "wf_human_snp"
    cpus 2
    input:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
    output:
        path "output/qual", emit: full_qual
    shell:
        '''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        mkdir output
        bgzip -fdc pileup.vcf.gz | \
        pypy $(which clair3.py) SelectQual \
                --output_fn output \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont 
        '''
}


process create_candidates {
    // Create BED files for candidate variants for "full alignment" network
    // from the previous full "pileup" variants across all chunks of all chroms
    //
    // Performed per chromosome; output a list of bed files one for each chunk.
    label "wf_human_snp"
    cpus 2
    input:
        each contig
        tuple path(ref), path(fai), path(ref_cache)
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectCandidates.py#L146
        path "candidate_bed/qual"
    output:
        tuple val(contig), path("candidate_bed/${contig}.*"), emit: candidate_bed, optional: true
    shell:
        // This creates BED files as candidate_bed/<ctg>.0_14 with candidates
        // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED
        // files.  All we really want are the BEDs, the file of filenames is
        // used for the purposes of parallel in the original workflow.
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218

        // TODO: would be nice to control the number of BEDs produced to enable
        // better parallelism.
        '''
        pypy $(which clair3.py) SelectCandidates \
            --pileup_vcf_fn pileup.vcf.gz \
            --split_folder candidate_bed \
            --ref_fn !{ref} \
            --var_pct_full !{params.var_pct_full} \
            --ref_pct_full !{params.ref_pct_full} \
            --platform ont \
            --ctgName !{contig}
        '''
}


process evaluate_candidates {
    // Run "full alignment" network for variants in a candidate bed file.
    // phased_bam just references the input BAM as it no longer contains phase information.
    label "wf_human_snp"
    cpus 1
    errorStrategy 'retry'
    input:
        tuple val(contig), path(phased_bam), path(phased_bam_index), path(phased_vcf)
        tuple val(contig), path(candidate_bed)
        tuple path(ref), path(fai), path(ref_cache)
        path(model)
    output:
        path "output/full_alignment_*.vcf", emit: full_alignment
    script:
        filename = candidate_bed.name
        def ref_path = "${ref_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
        """
        export REF_PATH=${ref_path}
        mkdir output
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn $model/full_alignment \
            --bam_fn $phased_bam \
            --call_fn output/full_alignment_${filename}.vcf \
            --sampleName ${params.sample_name} \
            --ref_fn ${ref} \
            --full_aln_regions ${candidate_bed} \
            --ctgName ${contig} \
            --add_indel_length \
            --gvcf ${params.GVCF} \
            --minMQ ${params.min_mq} \
            --minCoverage ${params.min_cov} \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --platform ont \
            --phased_vcf_fn ${phased_vcf}
        """
}


process aggregate_full_align_variants {
    // Sort and merge all "full alignment" variants
    label "wf_human_snp"
    cpus 2
    input:
        tuple path(ref), path(fai), path(ref_cache)
        path "full_alignment/*"
        path contigs
        path "gvcf_tmp_path/*"
    output:
        tuple path("full_alignment.vcf.gz"), path("full_alignment.vcf.gz.tbi"), emit: full_aln_vcf
        path "non_var.gvcf", optional: true, emit: non_var_gvcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir full_alignment \
            --output_fn full_alignment.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( bcftools index -n full_alignment.vcf.gz )" -eq 0 ]; then
            echo "[INFO] Exit in full-alignment variant calling"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ "!{params.GVCF}" == "true" ]; then
            pypy $(which clair3.py) SortVcf \
                --input_dir gvcf_tmp_path \
                --vcf_fn_suffix .tmp.gvcf \
                --output_fn non_var.gvcf \
                --sampleName !{params.sample_name} \
                --ref_fn !{ref} \
                --contigs_fn !{contigs}
        fi
        '''
}


process merge_pileup_and_full_vars{
    // Merge VCFs
    label "wf_human_snp"
    cpus 2
    input:
        each contig
        tuple path(ref), path(fai), path(ref_cache)
        tuple path(pile_up_vcf), path(pile_up_vcf_tbi)
        tuple path(full_aln_vcf), path(full_aln_vcf_tbi)
        path "non_var.gvcf"
        path "candidate_beds/*"
    output:
        tuple val(contig), path("output/merge_${contig}.vcf.gz"), path("output/merge_${contig}.vcf.gz.tbi"), emit: merged_vcf
        path "output/merge_${contig}.gvcf", optional: true, emit: merged_gvcf
    shell:
        '''
        mkdir output
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        pypy $(which clair3.py) MergeVcf \
            --pileup_vcf_fn !{pile_up_vcf} \
            --bed_fn_prefix candidate_beds \
            --full_alignment_vcf_fn !{full_aln_vcf} \
            --output_fn output/merge_!{contig}.vcf \
            --platform ont \
            --print_ref_calls False \
            --gvcf !{params.GVCF} \
            --haploid_precise False \
            --haploid_sensitive False \
            --gvcf_fn output/merge_!{contig}.gvcf \
            --non_var_gvcf_fn non_var.gvcf \
            --ref_fn !{ref} \
            --ctgName !{contig}

        bgzip -c output/merge_!{contig}.vcf > output/merge_!{contig}.vcf.gz
        tabix output/merge_!{contig}.vcf.gz
        '''
}


process post_clair_phase_contig {
    // Phase VCF for a contig
    label "wf_human_snp"
    cpus 4
    input:
        tuple val(contig), path(vcf), path(vcf_tbi), path(bam), path(bai), path(ref), path(fai), path(ref_cache)
    output:
        tuple val(contig), path("phased_${contig}.vcf.gz"), path("phased_${contig}.vcf.gz.tbi"), emit: vcf
    shell:
        '''
        if [[ "!{params.use_longphase}" == "true" ]]; then
            echo "Using longphase for phasing"
            # longphase needs decompressed 
            gzip -dc !{vcf} > variants.vcf
            longphase phase --ont -o phased_!{contig} \
                -s variants.vcf -b !{bam} -r !{ref} -t !{task.cpus}
            bgzip phased_!{contig}.vcf
        else
            echo "Using whatshap for phasing"
            whatshap phase \
                --output phased_!{contig}.vcf.gz \
                --reference !{ref} \
                --chromosome !{contig} \
                --distrust-genotypes \
                --ignore-read-groups \
                !{vcf} \
                !{bam}
        fi

        tabix -f -p vcf phased_!{contig}.vcf.gz

        '''
}


process aggregate_all_variants{
    label "wf_human_snp"
    cpus 4
    input:
        tuple path(ref), path(fai), path(ref_cache)
        path "merge_output/*"
        path "merge_outputs_gvcf/*"
        val(phase_vcf)
        path contigs
    output:
        tuple path("${params.sample_name}.wf_snp.vcf.gz"), path("${params.sample_name}.wf_snp.vcf.gz.tbi"), emit: final_vcf
        tuple path("${params.sample_name}.wf_snp.gvcf.gz"), path("${params.sample_name}.wf_snp.gvcf.gz.tbi"), emit: final_gvcf, optional: true
    shell:
        '''
        prefix="merge"
        phase_vcf=!{params.phase_vcf}
        if [[ $phase_vcf == "true" ]]; then
            prefix="phased"
        fi
        ls merge_output/*.vcf.gz | parallel --jobs 4 "bgzip -d {}"

        pypy $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --vcf_fn_prefix $prefix \
            --output_fn !{params.sample_name}.wf_snp.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( bgzip -fdc !{params.sample_name}.wf_snp.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
            echo "[INFO] Exit in all contigs variant merging"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ "!{params.GVCF}" == "true" ]; then
            pypy $(which clair3.py) SortVcf \
                --input_dir merge_outputs_gvcf \
                --vcf_fn_prefix merge \
                --vcf_fn_suffix .gvcf \
                --output_fn !{params.sample_name}.wf_snp.gvcf \
                --sampleName !{params.sample_name} \
                --ref_fn !{ref} \
                --contigs_fn !{contigs}
        fi

        echo "[INFO] Finish calling, output file: merge_output.vcf.gz"
        '''
}


process hap {
    label "happy"
    input:
        tuple path("clair.vcf.gz"), path("clair.vcf.gz.tbi")
        tuple path("ref.fasta"), path("ref.fasta.fai")
        path "truth.vcf"
        path "truth.bed"
    output:
        path "happy"
    shell:
        '''
        mkdir happy
        /opt/hap.py/bin/hap.py \
            truth.vcf \
            clair.vcf.gz \
            -f truth.bed \
            -r ref.fasta \
            -o happy \
            --engine=vcfeval \
            --threads=4 \
            --pass-only
        '''
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_snp {
    // publish inputs to output directory
    label "wf_human_snp"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


process getVersions {
    label "wf_human_snp"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        run_clair3.sh --version | sed 's/ /,/' >> versions.txt
        """
}


process getParams {
    label "wf_human_snp"
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


process vcfStats {
    label "wf_human_snp"
    cpus 1
    input:
        tuple path(vcf), path(index)
    output:
        file "variants.stats"
    """
    bcftools stats $vcf > variants.stats
    """
}


process makeReport {
    input:
        file read_summary
        tuple path(mosdepth_bed), path(mosdepth_dist), path(mosdepth_threshold) // MOSDEPTH_TUPLE
        file vcfstats
        path versions
        path "params.json"
    output:
        path "*report.html"
    script:
        report_name = "${params.sample_name}.wf-human-snp-report.html"
        wfversion = workflow.manifest.version
        if( workflow.commitId ){
            wfversion = workflow.commitId
        }
        """
        workflow-glue report \
        $report_name \
        --versions $versions \
        --params params.json \
        --read_stats $read_summary \
        --read_depth $mosdepth_dist \
        --vcf_stats $vcfstats \
        --revision $wfversion \
        --commit $workflow.commitId
        """
}


// This is a hilarious trick from CW to present models inside the container as
// outside the container, by exporting them out of the container back to workdir.
// This saves us passing around tuples of val(inside) and path(outside).
process lookup_clair3_model {
    label "wf_human_snp"
    input:
        path("lookup_table")
        val basecall_model
    output:
        path("model/")
    shell:
    '''
    clair3_model=$(resolve_clair3_model.py lookup_table '!{basecall_model}')
    cp -r ${CLAIR_MODELS_PATH}/${clair3_model} model
    echo "Basecall model: !{basecall_model}"
    echo "Clair3 model  : ${clair3_model}"
    '''
}
