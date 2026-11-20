#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd


PANEL_NAME = "1KG_Boehringer_ONT_Noyvert_888"

# Columnas UCSC de noyvert.bed, numeradas desde cero.
# Las coordenadas BED son 0-based, con extremo final exclusivo.
UCSC_COLUMNS = {
    "name": 3,
    "svtype": 9,
    "svlen": 10,
    "inslen": 11,
    "AF": 14,
    "AF_AFR": 15,
    "AF_AMR": 16,
    "AF_EAS": 17,
    "AF_EUR": 18,
    "AF_SAS": 19,
    "AF_UKB_IMPUTED": 27,
    "R2_LOO": 25,
    "CONCORDANCE_LOO": 26,
}

OUTPUT_COLUMNS = [
    "LR_DB",
    "LR_MATCH_STATUS",
    "LR_MATCH_ID",
    "LR_MATCH_TYPE",
    "LR_MATCH_COORD",
    "LR_MATCH_N_CANDIDATES",
    "LR_MATCH_RECIPROCAL_OVERLAP",
    "LR_MATCH_BP_DISTANCE",
    "LR_MATCH_LENGTH_DIFF_FRAC",
    "LR_AF",
    "LR_AF_AFR",
    "LR_AF_AMR",
    "LR_AF_EAS",
    "LR_AF_EUR",
    "LR_AF_SAS",
    "LR_AF_UKB_IMPUTED",
    "LR_R2_LOO",
    "LR_CONCORDANCE_LOO",
]


def normalize_chrom(value):
    chrom = str(value).strip()

    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    if chrom.upper() in ("M", "MT"):
        return "MT"

    return chrom.upper()


def parse_int(value):
    try:
        if value is None or str(value).strip() in ("", "NA", "."):
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def normalize_type(value):
    svtype = str(value).strip().upper()

    # Sniffles/AnnotSV pueden representar translocaciones como TRA;
    # UCSC registra el tipo correspondiente como BND.
    if svtype == "TRA":
        return "BND"

    return svtype


def load_reference(path):
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
    )

    if raw.shape[1] < 31:
        raise ValueError(
            "El BED tiene {} columnas; se esperaban al menos 31."
            .format(raw.shape[1])
        )

    ref = pd.DataFrame({
        "chrom": raw.iloc[:, 0].map(normalize_chrom),
        "start": pd.to_numeric(raw.iloc[:, 1], errors="raise"),
        "end": pd.to_numeric(raw.iloc[:, 2], errors="raise"),
        "name": raw.iloc[:, 3].astype(str),
        "type": raw.iloc[:, 9].map(normalize_type),
        "svlen": pd.to_numeric(raw.iloc[:, 10], errors="coerce"),
        "inslen": pd.to_numeric(raw.iloc[:, 11], errors="coerce"),
        "AF": raw.iloc[:, 14].astype(str),
        "AF_AFR": raw.iloc[:, 15].astype(str),
        "AF_AMR": raw.iloc[:, 16].astype(str),
        "AF_EAS": raw.iloc[:, 17].astype(str),
        "AF_EUR": raw.iloc[:, 18].astype(str),
        "AF_SAS": raw.iloc[:, 19].astype(str),
        "R2_LOO": raw.iloc[:, 25].astype(str),
        "CONCORDANCE_LOO": raw.iloc[:, 26].astype(str),
        "AF_UKB_IMPUTED": raw.iloc[:, 27].astype(str),
    })

    ref["start"] = ref["start"].astype(np.int64)
    ref["end"] = ref["end"].astype(np.int64)

    # Índice por cromosoma y tipo para no comparar cada SV con
    # las 107.445 variantes del panel.
    groups = {}

    for key, group in ref.groupby(["chrom", "type"], sort=False):
        groups[key] = {
            "start": group["start"].to_numpy(dtype=np.int64),
            "end": group["end"].to_numpy(dtype=np.int64),
            "svlen": group["svlen"].to_numpy(dtype=float),
            "inslen": group["inslen"].to_numpy(dtype=float),
            "name": group["name"].to_numpy(dtype=object),
            "AF": group["AF"].to_numpy(dtype=object),
            "AF_AFR": group["AF_AFR"].to_numpy(dtype=object),
            "AF_AMR": group["AF_AMR"].to_numpy(dtype=object),
            "AF_EAS": group["AF_EAS"].to_numpy(dtype=object),
            "AF_EUR": group["AF_EUR"].to_numpy(dtype=object),
            "AF_SAS": group["AF_SAS"].to_numpy(dtype=object),
            "AF_UKB_IMPUTED": group["AF_UKB_IMPUTED"].to_numpy(
                dtype=object
            ),
            "R2_LOO": group["R2_LOO"].to_numpy(dtype=object),
            "CONCORDANCE_LOO": group["CONCORDANCE_LOO"].to_numpy(
                dtype=object
            ),
        }

    chromosomes = set(ref["chrom"].unique())

    return groups, chromosomes, len(ref)


def empty_result(status):
    result = {
        "LR_DB": PANEL_NAME,
        "LR_MATCH_STATUS": status,
        "LR_MATCH_ID": "NA",
        "LR_MATCH_TYPE": "NA",
        "LR_MATCH_COORD": "NA",
        "LR_MATCH_N_CANDIDATES": "0",
        "LR_MATCH_RECIPROCAL_OVERLAP": "NA",
        "LR_MATCH_BP_DISTANCE": "NA",
        "LR_MATCH_LENGTH_DIFF_FRAC": "NA",
    }

    for col in [
        "LR_AF", "LR_AF_AFR", "LR_AF_AMR", "LR_AF_EAS",
        "LR_AF_EUR", "LR_AF_SAS", "LR_AF_UKB_IMPUTED",
        "LR_R2_LOO", "LR_CONCORDANCE_LOO",
    ]:
        result[col] = "NA"

    return result


def match_variant(row, groups, chromosomes, reciprocal_overlap,
                  insertion_distance, insertion_length_diff):
    chrom = normalize_chrom(row.get("SV_chrom", ""))
    svtype = normalize_type(row.get("SV_type", ""))

    start = parse_int(row.get("SV_start"))
    end = parse_int(row.get("SV_end"))
    length = parse_int(row.get("SV_length"))

    if start is None or end is None or not svtype:
        return empty_result("missing_coordinates")

    if chrom not in chromosomes:
        return empty_result("chromosome_not_in_panel")

    # Para BND/TRA se necesita comparar también el punto de ruptura
    # acompañante. No se fuerza una coincidencia basada en un solo punto.
    if svtype == "BND":
        return empty_result("BND_requires_mate_breakpoint")

    if svtype not in ("DEL", "DUP", "INV", "INS"):
        return empty_result("unsupported_type")

    key = (chrom, svtype)
    group = groups.get(key)

    if group is None:
        return empty_result("no_candidate_same_chrom_type")

    # AnnotSV informa la posición de inicio en coordenadas genómicas
    # 1-based. La normalizamos al inicio 0-based usado por UCSC.
    query_start = start - 1
    query_len = abs(length) if length is not None else abs(end - start)

    best_index = None
    candidate_count = 0
    result = empty_result("no_match_by_threshold")

    if svtype in ("DEL", "DUP", "INV"):
        if query_len <= 0:
            return empty_result("invalid_interval_length")

        query_end = query_start + query_len

        ref_start = group["start"]
        ref_end = group["end"]
        ref_len = ref_end - ref_start

        overlap = np.maximum(
            0,
            np.minimum(query_end, ref_end)
            - np.maximum(query_start, ref_start),
        )

        valid_ref = ref_len > 0
        fraction_query = overlap / float(query_len)

        fraction_ref = np.zeros(len(ref_len), dtype=float)
        fraction_ref[valid_ref] = (
            overlap[valid_ref] / ref_len[valid_ref]
        )

        scores = np.minimum(fraction_query, fraction_ref)
        valid = np.where(
            valid_ref & (scores >= reciprocal_overlap)
        )[0]

        candidate_count = len(valid)

        if candidate_count == 0:
            return empty_result("no_match_by_threshold")

        best_index = valid[np.argmax(scores[valid])]
        result = empty_result(
            "matched" if candidate_count == 1
            else "multiple_candidates_best_ranked"
        )

        result["LR_MATCH_RECIPROCAL_OVERLAP"] = (
            "{:.6f}".format(scores[best_index])
        )

    elif svtype == "INS":
        distances = np.abs(group["start"] - query_start)
        valid = distances <= insertion_distance

        ref_lengths = group["inslen"]
        size_diffs = np.full(len(ref_lengths), np.nan, dtype=float)

        if query_len > 0:
            known = np.isfinite(ref_lengths) & (ref_lengths > 0)
            size_diffs[known] = (
                np.abs(query_len - ref_lengths[known])
                / np.maximum(query_len, ref_lengths[known])
            )

            valid &= (
                ~known
                | (size_diffs <= insertion_length_diff)
            )

        candidate_indices = np.where(valid)[0]
        candidate_count = len(candidate_indices)

        if candidate_count == 0:
            return empty_result("no_match_by_threshold")

        # Prioriza proximidad y, cuando hay tamaño disponible,
        # similitud de longitud de la inserción.
        costs = (
            distances[candidate_indices] / float(insertion_distance)
        )

        known_diff = np.isfinite(size_diffs[candidate_indices])
        costs[known_diff] += size_diffs[candidate_indices][known_diff]

        best_index = candidate_indices[np.argmin(costs)]

        result = empty_result(
            "matched" if candidate_count == 1
            else "multiple_candidates_best_ranked"
        )

        result["LR_MATCH_BP_DISTANCE"] = str(
            int(distances[best_index])
        )

        if np.isfinite(size_diffs[best_index]):
            result["LR_MATCH_LENGTH_DIFF_FRAC"] = (
                "{:.6f}".format(size_diffs[best_index])
            )

    result["LR_MATCH_N_CANDIDATES"] = str(candidate_count)
    result["LR_MATCH_ID"] = str(group["name"][best_index])
    result["LR_MATCH_TYPE"] = svtype
    result["LR_MATCH_COORD"] = "{}:{}-{}".format(
        "chr" + chrom,
        int(group["start"][best_index]),
        int(group["end"][best_index]),
    )

    field_map = {
        "LR_AF": "AF",
        "LR_AF_AFR": "AF_AFR",
        "LR_AF_AMR": "AF_AMR",
        "LR_AF_EAS": "AF_EAS",
        "LR_AF_EUR": "AF_EUR",
        "LR_AF_SAS": "AF_SAS",
        "LR_AF_UKB_IMPUTED": "AF_UKB_IMPUTED",
        "LR_R2_LOO": "R2_LOO",
        "LR_CONCORDANCE_LOO": "CONCORDANCE_LOO",
    }

    for output_col, ref_col in field_map.items():
        value = group[ref_col][best_index]
        result[output_col] = (
            str(value) if str(value).strip() not in ("", ".", "NA")
            else "NA"
        )

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Añade frecuencias de UCSC Long-read SVs a AnnotSV."
    )
    parser.add_argument("--input", required=True, help="TSV original de AnnotSV")
    parser.add_argument("--reference", required=True, help="BED de UCSC noyvert.bb")
    parser.add_argument("--output", required=True, help="TSV enriquecido de salida")
    parser.add_argument("--reciprocal-overlap", type=float, default=0.70)
    parser.add_argument("--insertion-distance", type=int, default=100)
    parser.add_argument("--insertion-length-diff", type=float, default=0.50)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.is_file():
        raise SystemExit("No existe el TSV de entrada: {}".format(input_path))

    if not Path(args.reference).is_file():
        raise SystemExit("No existe el BED UCSC: {}".format(args.reference))

    print("Cargando panel UCSC...")
    groups, chromosomes, n_reference = load_reference(args.reference)
    print("Registros del panel:", n_reference)

    print("Leyendo TSV de AnnotSV...")
    table = pd.read_csv(
        input_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = ["SV_chrom", "SV_start", "SV_end", "SV_type"]
    missing = [col for col in required if col not in table.columns]

    if missing:
        raise SystemExit(
            "Faltan columnas en AnnotSV: {}".format(missing)
        )

    # Permite volver a ejecutar el script sobre un TSV enriquecido
    # sin duplicar las columnas LR_*.
    table = table.drop(
        [col for col in OUTPUT_COLUMNS if col in table.columns],
        axis=1,
        errors="ignore",
    )

    results = []
    cache = {}

    for row in table.to_dict("records"):
        key = (
            normalize_chrom(row.get("SV_chrom", "")),
            parse_int(row.get("SV_start")),
            parse_int(row.get("SV_end")),
            normalize_type(row.get("SV_type", "")),
            parse_int(row.get("SV_length")),
        )

        if key not in cache:
            cache[key] = match_variant(
                row,
                groups,
                chromosomes,
                args.reciprocal_overlap,
                args.insertion_distance,
                args.insertion_length_diff,
            )

        results.append(cache[key].copy())

    added = pd.DataFrame(results, index=table.index)
    final = pd.concat([table, added], axis=1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(output_path, sep="\t", index=False)

    counts = Counter(added["LR_MATCH_STATUS"])

    print("Salida:", output_path)
    print("Filas AnnotSV:", len(final))
    print("Estados de coincidencia:")
    for status, count in sorted(counts.items()):
        print("  {}: {}".format(status, count))


if __name__ == "__main__":
    main()
