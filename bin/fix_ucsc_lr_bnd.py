import argparse
import csv
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

# Permitir campos grandes, como la columna ANN de AnnotSV.
max_field_size = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_field_size)
        break
    except OverflowError:
        max_field_size //= 10

parser = argparse.ArgumentParser(
    description="Resuelve BND/TRA utilizando parejas candidatas UCSC."
)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--summary", type=Path, required=True)
parser.add_argument("--candidates", type=Path, required=True)
parser.add_argument("--reference", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--candidate-output", type=Path, required=True)
args = parser.parse_args()

src = args.input
summary_file = args.summary
candidates_file = args.candidates
bed_file = args.reference
out = args.output
candidate_out = args.candidate_output

for f in (src, summary_file, candidates_file, bed_file):
    if not f.is_file():
        raise SystemExit("Falta el fichero necesario: {}".format(f))


def norm_chrom(value):
    c = str(value).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "MT"
    return c.upper()


def pos_int(value):
    try:
        return str(int(float(value)))
    except (ValueError, TypeError):
        return str(value).strip()


def query_key(row):
    return (
        row.get("Samples_ID", ""),
        norm_chrom(row.get("SV_chrom", "")),
        pos_int(row.get("SV_start", "")),
        row.get("ALT", "")
    )


def parse_second_bp(text):
    m = re.search(r'([\[\]])([^:\[\]]+):(\d+)\1', str(text))
    if not m:
        return None
    return norm_chrom(m.group(2)), int(m.group(3))


def parse_primary_coord(text):
    m = re.match(r'^(?:chr)?([^:]+):(\d+)-(\d+)$', str(text))
    if not m:
        return None
    return norm_chrom(m.group(1)), int(m.group(2)), int(m.group(3))


def parse_mate_coord(text):
    m = re.match(r'^(?:chr)?([^:]+):(\d+)$', str(text))
    if not m:
        return None
    return norm_chrom(m.group(1)), int(m.group(2))


def read_tsv(path):
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames, rows


print("Leyendo informes BND existentes...")

summary_fields, summaries = read_tsv(summary_file)
candidate_fields, candidate_rows = read_tsv(candidates_file)

summary_map = {query_key(r): r for r in summaries}
candidates_by_query = defaultdict(list)

for r in candidate_rows:
    candidates_by_query[query_key(r)].append(r)

# Recuperar métricas adicionales desde el registro UCSC original.
reference_extra = {}

with open(bed_file, "r", newline="") as f:
    for line in f:
        b = line.rstrip("\n").split("\t")

        if len(b) < 28 or b[9].upper() != "BND":
            continue

        second = parse_second_bp(b[20])
        if second is None:
            continue

        ref_key = (
            b[3],
            norm_chrom(b[0]),
            int(b[1]),
            int(b[2]),
            second[0],
            second[1]
        )

        reference_extra[ref_key] = {
            "LR_AF_UKB_IMPUTED": b[27],
            "LR_R2_LOO": b[25],
            "LR_CONCORDANCE_LOO": b[26]
        }

print("Leyendo TSV enriquecido...")
input_fields, table = read_tsv(src)

output_fields = list(input_fields)

new_fields = [
    "LR_MATCH_MATE_COORD",
    "LR_MATCH_MATE_BP_DISTANCE",
    "LR_MATCH_ORIENTATION_MATCH"
]

af_fields = [
    "LR_AF",
    "LR_AF_AFR",
    "LR_AF_AMR",
    "LR_AF_EAS",
    "LR_AF_EUR",
    "LR_AF_SAS",
    "LR_AF_UKB_IMPUTED",
    "LR_R2_LOO",
    "LR_CONCORDANCE_LOO"
]

for c in new_fields + af_fields:
    if c not in output_fields:
        output_fields.append(c)


def empty_result(status):
    result = {
        "LR_MATCH_STATUS": status,
        "LR_MATCH_ID": "NA",
        "LR_MATCH_TYPE": "NA",
        "LR_MATCH_COORD": "NA",
        "LR_MATCH_MATE_COORD": "NA",
        "LR_MATCH_N_CANDIDATES": "0",
        "LR_MATCH_RECIPROCAL_OVERLAP": "NA",
        "LR_MATCH_BP_DISTANCE": "NA",
        "LR_MATCH_MATE_BP_DISTANCE": "NA",
        "LR_MATCH_ORIENTATION_MATCH": "NA",
        "LR_MATCH_LENGTH_DIFF_FRAC": "NA"
    }

    for c in af_fields:
        result[c] = "NA"

    return result


print("Procesando estados BND/TRA...")

unique_status = {}

for qkey, summary in summary_map.items():
    query_candidates = candidates_by_query.get(qkey, [])

    # Contar solo candidatas que coincidan también en orientación.
    compatible = [
        c for c in query_candidates
        if c.get("ORIENTATION_MATCH", "").strip().lower() == "true"
    ]

    original_status = summary.get("STATUS", "")

    status_map = {
        "one_candidate": "matched",
        "multiple_candidates": "BND_multiple_candidates",
        "orientation_conflict": "BND_orientation_conflict",
        "no_pair_within_distance": "BND_no_pair_within_distance",
        "ALT_unparseable": "BND_unparseable_ALT",
        "invalid_local_position": "BND_invalid_local_position"
    }

    status = status_map.get(
        original_status,
        "BND_" + original_status
    )

    # Asignar frecuencia únicamente cuando el diagnóstico inicial
    # y el recuento de candidatas compatibles concuerdan.
    chosen = None

    if original_status == "one_candidate" and len(compatible) == 1:
        chosen = compatible[0]
        status = "matched"
    elif original_status == "one_candidate":
        status = "BND_candidate_count_inconsistent"

    result = empty_result(status)
    result["LR_MATCH_N_CANDIDATES"] = str(len(compatible))

    if chosen is not None:
        primary = parse_primary_coord(
            chosen.get("UCSC_PRIMARY_COORD_BED", "")
        )
        mate = parse_mate_coord(
            chosen.get("UCSC_SECOND_BP", "")
        )

        if primary is None or mate is None:
            result = empty_result("BND_invalid_UCSC_coordinates")
        else:
            result.update({
                "LR_MATCH_ID": chosen.get("UCSC_ID", "NA"),
                "LR_MATCH_TYPE": "BND",
                "LR_MATCH_COORD": chosen.get(
                    "UCSC_PRIMARY_COORD_BED", "NA"
                ),
                "LR_MATCH_MATE_COORD": chosen.get(
                    "UCSC_SECOND_BP", "NA"
                ),
                "LR_MATCH_BP_DISTANCE": chosen.get(
                    "DIST_LOCAL_BP", "NA"
                ),
                "LR_MATCH_MATE_BP_DISTANCE": chosen.get(
                    "DIST_MATE_BP", "NA"
                ),
                "LR_MATCH_ORIENTATION_MATCH": "True",
                "LR_AF": chosen.get("LR_AF", "NA"),
                "LR_AF_AFR": chosen.get("LR_AF_AFR", "NA"),
                "LR_AF_AMR": chosen.get("LR_AF_AMR", "NA"),
                "LR_AF_EAS": chosen.get("LR_AF_EAS", "NA"),
                "LR_AF_EUR": chosen.get("LR_AF_EUR", "NA"),
                "LR_AF_SAS": chosen.get("LR_AF_SAS", "NA")
            })

            ref_key = (
                chosen.get("UCSC_ID", ""),
                primary[0],
                primary[1],
                primary[2],
                mate[0],
                mate[1]
            )

            extra = reference_extra.get(ref_key, {})

            for c in (
                "LR_AF_UKB_IMPUTED",
                "LR_R2_LOO",
                "LR_CONCORDANCE_LOO"
            ):
                result[c] = extra.get(c, "NA")

    unique_status[qkey] = result["LR_MATCH_STATUS"]

    # Actualizar todas las filas de AnnotSV pertenecientes a esta llamada.
    for row in table:
        if row.get("SV_type", "").strip().upper() not in ("TRA", "BND"):
            continue

        if query_key(row) != qkey:
            continue

        row.update(result)

# Guardar el TSV enriquecido nuevo.
with open(out, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=output_fields,
        delimiter="\t",
        extrasaction="ignore",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(table)

# Añadir al informe una columna que distingue la compatibilidad individual
# de cada candidata del estado global de la llamada.
final_candidate_fields = list(candidate_fields)

for c in (
    "CANDIDATE_ORIENTATION_COMPATIBLE",
    "CALL_SUMMARY_STATUS"
):
    if c not in final_candidate_fields:
        final_candidate_fields.append(c)

with open(candidate_out, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=final_candidate_fields,
        delimiter="\t",
        extrasaction="ignore",
        lineterminator="\n"
    )
    writer.writeheader()

    for r in candidate_rows:
        qkey = query_key(r)
        summary = summary_map.get(qkey, {})

        r["CANDIDATE_ORIENTATION_COMPATIBLE"] = (
            "True"
            if r.get("ORIENTATION_MATCH", "").strip().lower() == "true"
            else "False"
        )

        r["CALL_SUMMARY_STATUS"] = summary.get(
            "STATUS", "unknown"
        )

        writer.writerow(r)

print("\nTSV nuevo guardado:", out.resolve())
print("Informe de candidatas guardado:", candidate_out.resolve())

print("\nEstados por llamada BND/TRA única:")

counts = Counter(unique_status.values())

for status, count in sorted(counts.items()):
    print("  {}: {}".format(status, count))

print("\nCandidatas exportadas:", len(candidate_rows))
