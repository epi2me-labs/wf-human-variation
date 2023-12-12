
include { 
    haploblocks as haploblocks_joint
} from '../modules/local/common.nf'


process phase_all {
    // Phase VCF for a contig
    cpus 4
    input:
        tuple val(contig),
            path(xam),
            path(xam_idx),
            path(snp_vcf, stageAs: "SNP/snp.vcf.gz"),
            path(snp_vcf_tbi, stageAs: "SNP/snp.vcf.gz.tbi"),
            path(sv_vcf, stageAs: "SV/sv.vcf.gz"),
            path(sv_vcf_tbi, stageAs: "SV/sv.vcf.gz.tbi"),
            val(is_sv_vcf),
            path(ref),
            path(ref_idx),
            path(ref_cache),
            env(REF_PATH)
    output:
        path("phased_${contig}.vcf.gz"), emit: vcf
        path("phased_${contig}.vcf.gz.tbi"), emit: tbi
    script:
        // Define option for longphase
        def sv_option = is_sv_vcf ? "--sv-file sv.vcf" : ""
        def sv_preparation = is_sv_vcf ? "bcftools view -r ${contig} -O v ${sv_vcf} > sv.vcf" : ""
        // Longphase requires SNPs to perform joint phasing
        """
        echo "Using longphase for phasing"

        # Extract contig of interest
        bcftools view -r ${contig} -O v ${snp_vcf} > snp.vcf

        # Same for the SV if available, otherwise this line will be empty
        ${sv_preparation}

        # Run longphase
        longphase phase \
            --ont \
            --indels \
            -o tmp_${contig} \
            -s snp.vcf \
            -b ${xam} \
            -r ${ref} \
            -t ${task.cpus} ${sv_option}

        # Compress all new VCFs.
        bgzip tmp_${contig}.vcf && tabix -p vcf tmp_${contig}.vcf.gz
        if [ -e tmp_${contig}_SV.vcf ]; then
            bgzip tmp_${contig}_SV.vcf && tabix -p vcf tmp_${contig}_SV.vcf.gz
        fi

        # Combine with the SV, if present.
        if [ -e tmp_${contig}_SV.vcf.gz ]; then
            bcftools concat -a -O u tmp_${contig}.vcf.gz tmp_${contig}_SV.vcf.gz | bcftools sort -O z - > phased_${contig}.vcf.gz
            tabix -f -p vcf phased_${contig}.vcf.gz
        else
            mv tmp_${contig}.vcf.gz phased_${contig}.vcf.gz
            mv tmp_${contig}.vcf.gz.tbi phased_${contig}.vcf.gz.tbi
        fi
        """
}


process vcf_concat_all {
    // Phase VCF for a contig
    // Emit directly since it is the only output from this stage
    input:
        path(phased_vcfs, stageAs: "phased/*")
        path(phased_tbis, stageAs: "phased/*")
    output:
        tuple path("${params.sample_name}.wf_human_variation.phased.vcf.gz"), path("${params.sample_name}.wf_human_variation.phased.vcf.gz.tbi"), emit: combined
    script:
        """
        # Prepare correct input file
        bcftools concat -O u phased/*.vcf.gz | bcftools sort -O z - > ${params.sample_name}.wf_human_variation.phased.vcf.gz
        tabix -p vcf ${params.sample_name}.wf_human_variation.phased.vcf.gz
        """
}


workflow phasing {
    take:
        clair_vcf
        sv_vcf
        reference_ch
        haplotagged_bams
    main:
        // Prepare the input files to process
        // If SV requested, set the channel to execute with true
        if (params.sv){
            struct_vcf = sv_vcf.map{
                vcf, tbi -> [vcf, tbi, true]
            }
        // otherwise, instantiate a dummy channel set not to be processed with false
        } else {
            struct_vcf = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE").map{
                fname -> [fname, fname, false]
            }
        }

        // Phase everything
        haplotagged_bams
            .combine(clair_vcf)
            .combine(struct_vcf)
            .combine(reference_ch) | phase_all

        // Collect phased files
        vcfs = phase_all.out.vcf.collect()
        tbis = phase_all.out.tbi.collect()

        // Concatenate the phased files
        vcf_concat_all(vcfs, tbis)

        // Define joint haploblocks
        haploblocks_joint(vcf_concat_all.out, 'joint')

    emit:
        vcf_concat_all.out.combined.concat(haploblocks_joint.out.phase_blocks)
}
