process UCSC_LR_REFERENCE {

    tag 'UCSC hg38 Noyvert long-read SV reference'

    label 'process_medium'



    container "${params.ucsc_lr_container}"
    publishDir "${params.out_dir}/annot_SVs/with_ucsc_lr/reference",
        mode: 'copy'

    output:
    path 'noyvert.bb', emit: bigbed
    path 'noyvert.bed', emit: bed
    path 'noyvert.reference.sha256', emit: checksum

    script:
    """
    set -euo pipefail

    curl --fail --location --retry 3 --retry-all-errors \
        --connect-timeout 30 \
        "https://hgdownload.soe.ucsc.edu/gbdb/hg38/lrSv/noyvert.bb" \
        -o noyvert.bb

    test -s noyvert.bb

    bigBedToBed noyvert.bb noyvert.bed
    test -s noyvert.bed

    awk -F '\\t' '
        NF < 31 {
            printf "ERROR: fila %d con %d columnas; se necesitan al menos 31\\n", NR, NF > "/dev/stderr"
            bad = 1
            exit 1
        }
        END {
            if (NR == 0 || bad) exit 1
            printf "Filas BED validadas: %d\\n", NR
        }
    ' noyvert.bed

    sha256sum noyvert.bb noyvert.bed > noyvert.reference.sha256
    """
}

process ADD_UCSC_LR_FREQUENCY {

    tag "${annotated_tsv.name}"

    label 'process_medium'


    container "${params.ucsc_lr_container}"
    publishDir "${params.out_dir}",
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.SVs.annotated.UCSC_LR.tsv') ?
                "annot_SVs/with_ucsc_lr/${filename}" :
            filename.endsWith('.CNVs.annotated.UCSC_LR.tsv') ?
                "annot_CNVs/with_ucsc_lr/${filename}" :
            null
        }

    input:
    tuple val(meta), path(annotated_tsv)
    path reference_bed
    path matcher_script

    output:
    tuple val(meta),
        path("*.annotated.UCSC_LR.tsv"),
        emit: annotated_ucsc_lr

    script:
    def variant_type = meta.variant_source?.toString()

    if (!(variant_type in ['SVs', 'CNVs'])) {
        throw new IllegalArgumentException(
            "ADD_UCSC_LR_FREQUENCY: meta.variant_source debe ser 'SVs' o 'CNVs'. " +
            "Valor recibido: ${variant_type}"
        )
    }

    def suffix = ".${variant_type}.annotated.tsv"

    if (!annotated_tsv.name.endsWith(suffix)) {
        throw new IllegalArgumentException(
            "El TSV no termina en ${suffix}: ${annotated_tsv.name}"
        )
    }

    def sample_prefix = annotated_tsv.name.substring(
        0,
        annotated_tsv.name.length() - suffix.length()
    )

    def output_name =
        "${sample_prefix}.${variant_type}.annotated.UCSC_LR.tsv"

    """
    set -euo pipefail

    python "${matcher_script}" \
        --input "${annotated_tsv}" \
        --reference "${reference_bed}" \
        --output "${output_name}" \
        --reciprocal-overlap 0.70 \
        --insertion-distance 100 \
        --insertion-length-diff 0.50

    test -s "${output_name}"

    python - "${output_name}" <<'PY'
import csv
import sys

filename = sys.argv[1]
required = {
    "LR_MATCH_STATUS",
    "LR_MATCH_ID",
    "LR_AF",
    "LR_AF_EUR",
}

with open(filename, "r", encoding="utf-8", newline="") as handle:
    reader = csv.reader(handle, delimiter="\\t")
    header = next(reader, None)

if not header:
    raise SystemExit("ERROR: el TSV final está vacío")

missing = sorted(required.difference(header))
if missing:
    raise SystemExit(
        "ERROR: faltan columnas UCSC: " + ", ".join(missing)
    )

print("Validación de cabecera correcta:", filename)
PY
    """
}