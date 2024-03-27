import groovy.json.JsonBuilder

def phaser_memory = params.use_longphase ? [8.GB, 32.GB, 56.GB] : [4.GB, 8.GB, 12.GB]
def haptag_memory = [4.GB, 8.GB, 12.GB]

// As of Clair3 v1.0.6, set `--min_snp_af` and `--min_indel_af` to 0 with `--vcf_fn`.
def snp_min_af = params.vcf_fn ? "--snp_min_af 0.0": "--snp_min_af ${params.snp_min_af}"
def indel_min_af = params.vcf_fn ? "--indel_min_af 0.0" : "--indel_min_af ${params.indel_min_af}"

process make_chunks {
    // Do some preliminaries. Ordinarily this would setup a working directory
    // that all other commands would make use off, but all we need here are the
    // list of contigs and chunks.
    label "wf_human_snp"
    cpus 1
    memory 4.GB
    input:
        tuple path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path bed
        path model_path
        val chromosome_codes
        path genotyping_vcf, stageAs: "genotyping_vcf/"
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
        path "clair_output/tmp/CMD", emit: cmd_file
        path "clair_output/tmp/split_beds", emit: split_beds, optional: true
    script:
        // Prepare input BED and genotyping VCF arguments (mutually exclusive, checked in workflow)
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed}" : ''
        def bedprnt = bed.name != 'OPTIONAL_FILE' ? "--bed_fn=${bed}" : ''
        def vcfargs = genotyping_vcf.baseName != "OPTIONAL_FILE" ? "--vcf_fn ${genotyping_vcf}" : ""
        def vcfprnt = genotyping_vcf.baseName != "OPTIONAL_FILE" ? "--vcf_fn=${genotyping_vcf}" : ""
        // Define contigs in order to enforce the mitochondrial genome calling, which is otherwise skipped.
        String ctgs = chromosome_codes.join(',')
        def ctg_name = "--ctg_name ${ctgs}"
        // If a single contig is required, then set it as option
        if (params.ctg_name){
            ctg_name = "--ctg_name ${params.ctg_name}"
        }
        // If all contigs are required, then set the ctg_name to EMPTY
        if (params.include_all_ctgs || bed.name != 'OPTIONAL_FILE'){
            ctg_name = '--ctg_name "EMPTY"'
        }
        """
        # CW-2456: save command line to add to VCF file (very long command...)
        mkdir -p clair_output/tmp
        echo "run_clair3.sh --bam_fn=${xam} ${bedprnt} --ref_fn=${ref} ${vcfprnt} --output=clair_output --platform=ont --sample_name=${params.sample_name} --model_path=${model_path.simpleName} --ctg_name=${params.ctg_name} ${ctg_name} --include_all_ctgs=${params.include_all_ctgs} --chunk_num=0 --chunk_size=5000000 --qual=${params.min_qual} --var_pct_full=${params.var_pct_full} --ref_pct_full=${params.ref_pct_full} ${snp_min_af} ${indel_min_af} --min_contig_size=${params.min_contig_size}" > clair_output/tmp/CMD
        # CW-2456: prepare other inputs normally
        python \$(which clair3.py) CheckEnvs \
            --bam_fn ${xam} \
            ${bedargs} \
            --output_fn_prefix clair_output \
            --ref_fn ${ref} \
            ${vcfargs} \
            ${ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs ${params.include_all_ctgs} \
            --threads 1  \
            --qual ${params.min_qual} \
            --sampleName ${params.sample_name} \
            --var_pct_full ${params.var_pct_full} \
            --ref_pct_full ${params.ref_pct_full} \
            ${snp_min_af} \
            ${indel_min_af} \
            --min_contig_size ${params.min_contig_size} \
            --cmd_fn clair_output/tmp/CMD
        """
}


process pileup_variants {
    // Calls variants per region ("chunk") using pileup network.
    label "wf_human_snp"
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 1
    input:
        each region
        tuple path(xam), path(xam_idx)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path model
        path bed
        path command
        path "split_bed"
    output:
        // TODO: make this explicit, why is pileup VCF optional?
        path "pileup_*.vcf", optional: true, emit: pileup_vcf_chunks
        path "gvcf_tmp_path/*", optional: true, emit: pileup_gvcf_chunks
    script:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        // note: snp_min_af and indel_min_af have an impact on performance
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed} --extend_bed split_bed/${region.contig}" : ''
        """
        python \$(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn ${model}/pileup \
            --bam_fn ${xam} \
            --call_fn pileup_${region.contig}_${region.chunk_id}.vcf \
            --ref_fn ${ref} \
            --ctgName ${region.contig} \
            --chunk_id ${region.chunk_id} \
            --chunk_num ${region.total_chunks} \
            --platform ont \
            --fast_mode False \
            ${snp_min_af} \
            ${indel_min_af} \
            --minMQ ${params.min_mq} \
            --minCoverage ${params.min_cov} \
            --call_snp_only False \
            --gvcf ${params.GVCF} \
            --base_err ${params.base_err} \
            --gq_bin_size ${params.gq_bin_size} \
            --temp_file_dir gvcf_tmp_path \
            --cmd_fn ${command} \
            --pileup \
            ${bedargs}
        """
}


process aggregate_pileup_variants {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing.
    label "wf_human_snp"
    cpus 2
    memory 4.GB
    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        // these need to be named as original, as program uses info from
        // contigs file to filter
        path "input_vcfs/*"
        path contigs
        path command
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
            --contigs_fn !{contigs} \
            --cmd_fn !{command}

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
    memory 4.GB
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

process phase_contig {
    // Tags reads in an input BAM from heterozygous SNPs
    // The haplotag step was removed in clair-v0.1.11 so this step re-emits
    //   the original BAM and BAI as phased_bam for compatability,
    //   but adds the VCF as it is now tagged with phasing information
    //   used later in the full-alignment model
    cpus 4
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple val(contig), path(het_snps), path(het_snps_tbi), path(xam), path(xam_idx), path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(contig), path(xam), path(xam_idx), path("phased_${contig}.vcf.gz"), emit: phased_bam_and_vcf
    script:
        if (params.use_longphase)
        """
        echo "Using longphase for phasing"
        longphase phase --ont -o phased_${contig} \
            -s ${het_snps} -b ${xam} -r ${ref} -t ${task.cpus}
        bgzip phased_${contig}.vcf
        tabix -f -p vcf phased_${contig}.vcf.gz
        """
        else
        """
        echo "Using whatshap for phasing"
        whatshap phase \
            --output phased_${contig}.vcf.gz \
            --reference ${ref} \
            --chromosome ${contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            ${het_snps} \
            ${xam}

        tabix -f -p vcf phased_${contig}.vcf.gz
        """
}

process cat_haplotagged_contigs {
    label "wf_human_snp"
    cpus 4
    memory 15.GB // cat should not need this, but weirdness occasionally strikes
    input:
        path contig_bams // intermediate input always BAM here
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        tuple val(xam_fmt), val(xai_fmt)
    output:
        tuple path("${params.sample_name}.haplotagged.${xam_fmt}"), path("${params.sample_name}.haplotagged.${xam_fmt}.${xai_fmt}"), emit: merged_xam
    script:
    def threads = Math.max(task.cpus - 1, 1)
    """
    # ensure this bit is idempotent as it will inevitably not be called so
    if [ -f seq_list.txt ]; then
        rm seq_list.txt
    fi
    # pick the "first" bam and read its SQ list to determine sort order
    samtools view -H --no-PG `ls *_hp.bam | head -n1` | grep '^@SQ' | sed -nE 's,.*SN:([^[:space:]]*).*,\\1,p' > seq_list.txt
    # append present contigs to a file of file names, to cat according to SQ order
    while read sq; do
        if [ -f "\${sq}_hp.bam" ]; then
            echo "\${sq}_hp.bam" >> cat.fofn
        fi
    done < seq_list.txt
    if [ ! -s cat.fofn ]; then
        echo "No haplotagged inputs to cat? Are the input file names staged correctly?"
        exit 70 # EX_SOFTWARE
    fi

    # cat just cats, if we want bam, we'll have to deal with that ourselves
    if [ "${xam_fmt}" = "cram" ]; then
        samtools cat -b cat.fofn --no-PG -o - | samtools view --no-PG -@ ${threads} --reference ${ref} -O CRAM --write-index -o "${params.sample_name}.haplotagged.cram##idx##${params.sample_name}.haplotagged.cram.crai"
    else
        samtools cat -b cat.fofn --no-PG -@ ${threads} -o "${params.sample_name}.haplotagged.bam"
        samtools index -@ ${threads} -b "${params.sample_name}.haplotagged.bam"
    fi
    """
}

process get_qual_filter {
    // Determines quality filter for selecting candidate variants for second
    // stage "full alignment" calling.
    label "wf_human_snp"
    cpus 2
    memory 4.GB
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
    memory 4.GB
    input:
        each contig
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
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
    memory { 8.GB * task.attempt }
    maxRetries 2
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple val(contig), path(phased_xam), path(phased_xam_idx), path(phased_vcf)
        tuple val(contig), path(candidate_bed)
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path(model)
        path(command)
    output:
        path "output/full_alignment_*.vcf", emit: full_alignment
    script:
        filename = candidate_bed.name
        """
        mkdir output
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVariantsFromCffi \
            --chkpnt_fn ${model}/full_alignment \
            --bam_fn ${phased_xam} \
            --call_fn output/full_alignment_${filename}.vcf \
            --sampleName ${params.sample_name} \
            --ref_fn ${ref} \
            --full_aln_regions ${candidate_bed} \
            --ctgName ${contig} \
            --add_indel_length \
            --gvcf ${params.GVCF} \
            --minMQ ${params.min_mq} \
            --minCoverage ${params.min_cov} \
            ${snp_min_af} \
            ${indel_min_af} \
            --platform ont \
            --cmd_fn ${command} \
            --phased_vcf_fn ${phased_vcf}
        """
}


process aggregate_full_align_variants {
    // Sort and merge all "full alignment" variants
    label "wf_human_snp"
    cpus 2
    memory { 4.GB * task.attempt }
    maxRetries 2
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path "full_alignment/*"
        path contigs
        path "gvcf_tmp_path/*"
        path command
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
            --cmd_fn !{command} \
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
                --cmd_fn !{command} \
                --contigs_fn !{contigs}
        fi
        '''
}


process merge_pileup_and_full_vars{
    // Merge VCFs
    label "wf_human_snp"
    cpus 2
    memory 4.GB
    input:
        each contig
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
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
    // CW-2383: now uses base image to allow phasing of both snps and indels
    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple val(contig), 
            path(vcf), path(vcf_tbi), 
            path(xam), path(xam_idx), 
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(contig), 
            path("phased_${contig}.vcf.gz"), path("phased_${contig}.vcf.gz.tbi"), 
            emit: vcf
        tuple val(contig), 
            path("phased_${contig}.vcf.gz"), path("phased_${contig}.vcf.gz.tbi"), 
            path(xam), path(xam_idx), 
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH),
            emit: for_tagging
    script:
    if (params.use_longphase)
        """
        echo "Using longphase for phasing"
        longphase phase --ont -o phased_${contig} \
            -s ${vcf} -b ${xam} -r ${ref} -t ${task.cpus}
        bgzip phased_${contig}.vcf
        tabix -f -p vcf phased_${contig}.vcf.gz
    """
    else
    """
        # REF_PATH points to the reference cache and allows faster parsing of CRAM files
        echo "Using whatshap for phasing"
        whatshap phase \
            --output phased_${contig}.vcf.gz \
            --reference ${ref} \
            --chromosome ${contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            ${vcf} \
            ${xam}
        tabix -f -p vcf phased_${contig}.vcf.gz
        """
}

process post_clair_contig_haplotag {
    // Tags reads in an input BAM from heterozygous SNPs
    // Also haplotag for those modes that need it
    // We emit BAM as the STR workflow does not fully support CRAM, and so the
    // STR workflow can start while the haplotagged XAM is being catted and
    // written for the final output

    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { haptag_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple val(contig),
            path(vcf), path(tbi),
            path(xam), path(xam_idx),
            path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
    output:
        tuple val(contig), path("${contig}_hp.bam"), path("${contig}_hp.bam.bai"), emit: phased_bam
    script:
    """
    whatshap haplotag \
        --reference ${ref} \
        --ignore-read-groups \
        --regions ${contig} \
        phased_${contig}.vcf.gz \
        ${xam} \
    | samtools view -O bam --reference $ref -@3 -o ${contig}_hp.bam##idx##${contig}_hp.bam.bai --write-index
    """
}


process aggregate_all_variants{
    label "wf_human_snp"
    cpus 4
    memory { 8.GB * task.attempt }
    maxRetries 2
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH)
        path "merge_output/*"
        path "merge_outputs_gvcf/*"
        val(phase_vcf)
        path contigs
        path command
    output:
        tuple path("${params.sample_name}.wf_snp.vcf.gz"), path("${params.sample_name}.wf_snp.vcf.gz.tbi"), emit: final_vcf
        tuple path("${params.sample_name}.wf_snp.gvcf.gz"), path("${params.sample_name}.wf_snp.gvcf.gz.tbi"), emit: final_gvcf, optional: true
    script:
        def prefix = params.phased || params.str ? "phased" : "merge"
        """
        ls merge_output/*.vcf.gz | parallel --jobs 4 "bgzip -d {}"

        pypy \$(which clair3.py) SortVcf \
            --input_dir merge_output \
            --vcf_fn_prefix $prefix \
            --output_fn ${params.sample_name}.wf_snp.vcf \
            --sampleName ${params.sample_name} \
            --ref_fn ${ref} \
            --cmd_fn ${command} \
            --contigs_fn ${contigs}

        if [ "\$( bgzip -fdc ${params.sample_name}.wf_snp.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
            echo "[INFO] Exit in all contigs variant merging"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ "${params.GVCF}" == "true" ]; then
            pypy \$(which clair3.py) SortVcf \
                --input_dir merge_outputs_gvcf \
                --vcf_fn_prefix merge \
                --vcf_fn_suffix .gvcf \
                --output_fn tmp.gvcf \
                --sampleName ${params.sample_name} \
                --ref_fn ${ref} \
                --cmd_fn ${command} \
                --contigs_fn ${contigs}

                # Reheading samples named "SAMPLE" to params.sample_name. If no 
                echo "SAMPLE" "${params.sample_name}" > rename.txt
                bcftools reheader -s rename.txt tmp.gvcf.gz > ${params.sample_name}.wf_snp.gvcf.gz
                bcftools index -t ${params.sample_name}.wf_snp.gvcf.gz && rm tmp.gvcf.gz rename.txt
        fi

        echo "[INFO] Finish calling, output file: merge_output.vcf.gz"
        """
}


process refine_with_sv {
    label "wf_human_snp"
    cpus 4
    memory { 8.GB * task.attempt - 1.GB }
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}

    input:
        tuple path(ref), path(ref_idx), path(ref_cache), env(REF_PATH) 
        tuple path(clair_vcf, stageAs: 'clair.vcf.gz'), path(clair_tbi, stageAs: 'clair.vcf.gz.tbi'), val(contig)
        tuple path(xam), path(xam_idx), val(meta) // this may be a haplotagged_bam or input CRAM 
        path sniffles_vcf
    output:
        tuple path("${params.sample_name}.${contig}.wf_snp.vcf.gz"), path("${params.sample_name}.${contig}.wf_snp.vcf.gz.tbi"), emit: final_vcf
    shell:
        '''
        pypy $(which clair3.py) SwitchZygosityBasedOnSVCalls \\
            --bam_fn !{xam} \\
            --clair3_vcf_input clair.vcf.gz \\
            --sv_vcf_input !{sniffles_vcf} \\
            --vcf_output '!{params.sample_name}.!{contig}.wf_snp.vcf' \\
            --threads !{task.cpus} \\
            --ctg_name '!{contig}'
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
    cpus 2
    input:
        tuple path(vcf), path(index)
    output:
        file "variants.stats"
    """
    bcftools stats --threads ${task.cpus - 1} $vcf > variants.stats
    """
}


process makeReport {
    cpus 1
    memory 16.GB
    input:
        file vcfstats
        path versions
        path "params.json"
        path clinvar_vcf
    output:
        path "*report.html"
    script:
        def clinvar = clinvar_vcf ?: ""
        def annotation = params.annotation ? "" : "--skip_annotation"

        report_name = "${params.sample_name}.wf-human-snp-report.html"
        wfversion = workflow.manifest.version
        if( workflow.commitId ){
            wfversion = workflow.commitId
        }
        """
        workflow-glue report_snp \
        $report_name \
        --versions $versions \
        --params params.json \
        --vcf_stats $vcfstats \
        --sample_name $params.sample_name \
        --clinvar_vcf $clinvar \
        $annotation
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
        tuple env(clair3_model), path("model/")
    shell:
    '''
    clair3_model=$(resolve_clair3_model.py lookup_table '!{basecall_model}')
    cp -r ${CLAIR_MODELS_PATH}/${clair3_model} model
    echo "Basecall model: !{basecall_model}"
    echo "Clair3 model  : ${clair3_model}"
    '''
}
