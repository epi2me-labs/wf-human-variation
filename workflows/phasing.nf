
include { 
    haploblocks as haploblocks_joint;
    concat_vcfs as concat_phased_vcfs;
} from '../modules/local/common.nf'

def phaser_memory = [8.GB, 32.GB, 56.GB]

process phase_all {
    // Phase VCF for a contig
    cpus 4
    // Define memory from phasing tool and number of attempt
    memory { phaser_memory[task.attempt - 1] }
    maxRetries 2
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
        """
        echo "Using longphase for phasing"

        # Extract contig of interest
        bcftools view --threads ${task.cpus} -r ${contig} -O v ${snp_vcf} > snp.vcf

        # Same for the SV if available, otherwise this line will be empty
        bcftools view --threads ${task.cpus} -r ${contig} -O v ${sv_vcf} > sv.vcf

        # Run longphase
        longphase phase \
            --ont \
            -o tmp_${contig} \
            -s snp.vcf \
            --sv-file sv.vcf \
            -b ${xam} \
            -r ${ref} \
            -t ${task.cpus}

        # Longphase adds a longphaseVersion header whenever a `#` is found in the input file
        # See https://github.com/twolinin/longphase/issues/30
        # Remove these spurious lines, then compress and index the files.
        awk 'BEGIN{isheader=1}; \$0~"##" && isheader==1 {print}; \$1=="#CHROM" {isheader=0; print}; \$1!~"#" {print}' tmp_${contig}.vcf \
        | bgzip -c > tmp_${contig}.vcf.gz && tabix -p vcf tmp_${contig}.vcf.gz

        # Repeat the same, but for the SVs.
        awk 'BEGIN{isheader=1}; \$0~"##" && isheader==1 {print}; \$1=="#CHROM" {isheader=0; print}; \$1!~"#" {print}' tmp_${contig}_SV.vcf \
        | bgzip -c > tmp_${contig}_SV.vcf.gz && tabix -p vcf tmp_${contig}_SV.vcf.gz
        
        # Combine with the SV.
        bcftools concat --threads ${task.cpus} -a -O u tmp_${contig}.vcf.gz tmp_${contig}_SV.vcf.gz | bcftools sort -O z - > phased_${contig}.vcf.gz
        tabix -f -p vcf phased_${contig}.vcf.gz
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
