
process mahi_mahi {
    tag { meta.alias }

    container "${params.mahi_mahi_container}"

    cpus 4
    memory 16.GB

    publishDir "${params.out_dir}/mahi-mahi",
        mode: "copy",
        overwrite: true

    input:
    tuple path(xam), path(xai), val(meta)
    tuple path(ref), path(ref_idx)
    val(sex)

    output:
    tuple path("${meta.alias}.vcf"),
          path("${meta.alias}.tsv"),
          path("${meta.alias}.tagged.bam"),
          path("${meta.alias}.config.yaml"),
          path("${meta.alias}.log"),
          emit: results

    path("${meta.alias}_plots",
         type: 'dir',
         optional: true,
         emit: plots)

    script:
    def sexValue = sex?.toString()?.trim()

    if (!(sexValue in ['XX', 'XY'])) {
        error "Mahi-mahi requires sex XX or XY; received: '${sexValue}'"
    }

    def xamFlag = xam.name.toLowerCase().endsWith('.cram') \
        ? '--cram' \
        : '--bam'

    def repeatId = params.mahi_mahi_repeat_id?.toString()?.trim()

    if (repeatId && !(repeatId ==~ /[A-Za-z0-9_.-]+/)) {
        error "Invalid Mahi-mahi repeat ID: '${repeatId}'"
    }

    def repeatArg = repeatId ? "--repeat-id ${repeatId}" : ""

    """
    set -euo pipefail

    echo "Mahi-mahi version:"
    mahi-mahi --version

    mahi-mahi \
        ${xamFlag} "${xam}" \
        --ref "${ref}" \
        --loci hg38 \
        --sample "${meta.alias}" \
        --sex "${sexValue}" \
        --out-dir mahi_out \
        ${repeatArg}

    shopt -s nullglob

for f in mahi_out/*; do
    mv -- "\$f" .
done

rmdir mahi_out
    """
}
