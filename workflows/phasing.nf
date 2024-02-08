
include { 
    haploblocks as haploblocks_joint;
    concat_vcfs as concat_phased_vcfs;
} from '../modules/local/common.nf'

def phaser_memory = params.use_longphase ? [8.GB, 32.GB, 56.GB] : [4.GB, 8.GB, 12.GB]

process phase_all {
    // Phase VCF for a contig
    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 3
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
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
        tuple path("phased_${contig}.vcf.gz"), path("phased_${contig}.vcf.gz.tbi"), emit: phased_vcf
    script:
        // Define option for longphase
        def sv_option = is_sv_vcf ? "--sv-file sv.vcf" : ""
        def sv_preparation = is_sv_vcf ? "bcftools view --threads ${task.cpus} -r ${contig} -O v ${sv_vcf} > sv.vcf" : ""
        // Longphase requires SNPs to perform joint phasing
        """
        echo "Using longphase for phasing"

        # Extract contig of interest
        bcftools view --threads ${task.cpus} -r ${contig} -O v ${snp_vcf} > snp.vcf

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
            bcftools concat --threads ${task.cpus} -a -O u tmp_${contig}.vcf.gz tmp_${contig}_SV.vcf.gz | bcftools sort -O z - > phased_${contig}.vcf.gz
            tabix -f -p vcf phased_${contig}.vcf.gz
        else
            mv tmp_${contig}.vcf.gz phased_${contig}.vcf.gz
            mv tmp_${contig}.vcf.gz.tbi phased_${contig}.vcf.gz.tbi
        fi
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

        // Concatenate the phased files
        phased_vcf = concat_phased_vcfs(
            phase_all.out.phased_vcf.collect(),
            "${params.sample_name}.wf_human_variation.phased"
        ).final_vcf

        // Define joint haploblocks
        haploblocks_joint(phased_vcf, 'human_variation')

    emit:
        phased_vcf.concat(haploblocks_joint.out.phase_blocks)
}
