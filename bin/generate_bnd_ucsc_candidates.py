#!/usr/bin/env python3

import argparse
import csv
import re


def normalize_second_bp(value):
    """Devuelve únicamente cromosoma:posición del segundo breakend."""
    match = re.search(r"([^:\[\]]+):(\d+)", value or "")
    if match:
        return f"{match.group(1)}:{match.group(2)}"
    return value
import sys
from pathlib import Path
from collections import Counter, defaultdict

LIMIT = sys.maxsize
while True:
    try:
        csv.field_size_limit(LIMIT)
        break
    except OverflowError:
        LIMIT //= 10

BREAKEND_RE = re.compile(r'([\[\]])([^:\[\]]+):(\d+)\1')


def norm_chrom(value):
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in ("M", "MT"):
        return "MT"
    return chrom.upper()


def pos_int(value):
    try:
        return str(int(float(value)))
    except (ValueError, TypeError):
        return str(value).strip()


def parse_breakend(text):
    match = BREAKEND_RE.search(str(text))
    if not match:
        return None

    # VCF/BED ALT-style notation:
    # sequence before the brackets -> local orientation '+'
    # sequence after the brackets  -> local orientation '-'
    # '[' means remote '+', ']' means remote '-'
    local_orientation = "+" if match.start() > 0 else "-"
    remote_orientation = "+" if match.group(1) == "[" else "-"

    return (
        norm_chrom(match.group(2)),
        int(match.group(3)),
        (local_orientation, remote_orientation),
    )


def query_key(row):
    return (
        row.get("Samples_ID", ""),
        norm_chrom(row.get("SV_chrom", "")),
        pos_int(row.get("SV_start", "")),
        row.get("ALT", ""),
    )


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruye las candidatas UCSC para llamadas BND/TRA."
    )
    parser.add_argument("--input", required=True, help="TSV AnnotSV enriquecido con UCSC")
    parser.add_argument("--reference", required=True, help="BED UCSC noyvert")
    parser.add_argument("--summary", required=True, help="TSV de resumen BND")
    parser.add_argument("--candidates", required=True, help="TSV de parejas candidatas")
    parser.add_argument("--max-distance", type=int, default=100)
    args = parser.parse_args()

    source = Path(args.input)
    reference_path = Path(args.reference)

    if not source.is_file():
        raise SystemExit(f"No existe el TSV de entrada: {source}")
    if not reference_path.is_file():
        raise SystemExit(f"No existe la referencia UCSC: {reference_path}")
    if args.max_distance < 0:
        raise SystemExit("--max-distance no puede ser negativo")

    # Leer una sola consulta por llamada única BND/TRA.
    queries = {}
    required = {"Samples_ID", "SV_chrom", "SV_start", "SV_type", "ALT"}

    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Faltan columnas de entrada: {sorted(missing)}")

        for row in reader:
            if row.get("SV_type", "").strip().upper() not in ("BND", "TRA"):
                continue
            queries.setdefault(query_key(row), row)

    # Indexar las variantes UCSC por sus dos cromosomas.
    reference_index = defaultdict(list)
    n_reference = 0

    with reference_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 31 or fields[9].strip().upper() != "BND":
                continue

            parsed_mate = parse_breakend(fields[20])
            if parsed_mate is None:
                continue

            try:
                start0 = int(fields[1])
                end0 = int(fields[2])
            except ValueError:
                continue

            record = {
                "chrom_raw": fields[0],
                "chrom": norm_chrom(fields[0]),
                "start0": start0,
                "end0": end0,
                # Convertir BED 0-based a posición de ruptura 1-based.
                "primary_pos": start0 + 1,
                "name": fields[3],
                "second_bp": fields[20],
                "mate_chrom": parsed_mate[0],
                "mate_pos": parsed_mate[1],
                "orientation": parsed_mate[2],
                "AF": fields[14],
                "AF_AFR": fields[15],
                "AF_AMR": fields[16],
                "AF_EAS": fields[17],
                "AF_EUR": fields[18],
                "AF_SAS": fields[19],
            }

            reference_index[
                (record["chrom"], record["mate_chrom"])
            ].append(record)
            n_reference += 1

    summary_fields = [
        "Samples_ID", "SV_chrom", "SV_start", "SV_type", "ALT",
        "MATE_CHROM", "MATE_POS", "STATUS",
        "N_POSITION_CANDIDATES", "N_ORIENTATION_MATCHES",
    ]

    candidate_fields = [
        "Samples_ID", "SV_chrom", "SV_start", "SV_type", "ALT",
        "SAMPLE_MATE_CHROM", "SAMPLE_MATE_POS",
        "UCSC_ID", "UCSC_PRIMARY_COORD_BED", "UCSC_SECOND_BP",
        "ENDPOINT_ORDER", "DIST_LOCAL_BP", "DIST_MATE_BP",
        "ORIENTATION_MATCH",
        "LR_AF", "LR_AF_AFR", "LR_AF_AMR",
        "LR_AF_EAS", "LR_AF_EUR", "LR_AF_SAS", "STATUS",
    ]

    summaries = []
    candidate_rows = []

    # Orden determinista para poder comparar con los informes anteriores.
    for qkey in sorted(queries):
        row = queries[qkey]
        raw_chrom = row.get("SV_chrom", "")
        raw_start = row.get("SV_start", "")
        svtype = row.get("SV_type", "").strip().upper()
        alt = row.get("ALT", "")
        sample = row.get("Samples_ID", "")

        parsed_alt = parse_breakend(alt)

        base = {
            "Samples_ID": sample,
            "SV_chrom": raw_chrom,
            "SV_start": raw_start,
            "SV_type": svtype,
            "ALT": alt,
            "MATE_CHROM": "NA",
            "MATE_POS": "NA",
        }

        if parsed_alt is None:
            summaries.append({
                **base,
                "STATUS": "ALT_unparseable",
                "N_POSITION_CANDIDATES": "0",
                "N_ORIENTATION_MATCHES": "0",
            })
            continue

        mate_chrom, mate_pos, query_orientation = parsed_alt
        base["MATE_CHROM"] = mate_chrom
        base["MATE_POS"] = str(mate_pos)

        try:
            local_pos = int(float(raw_start))
        except (TypeError, ValueError):
            summaries.append({
                **base,
                "STATUS": "invalid_local_position",
                "N_POSITION_CANDIDATES": "0",
                "N_ORIENTATION_MATCHES": "0",
            })
            continue

        local_chrom = norm_chrom(raw_chrom)
        position_candidates = []

        # Direct: local del paciente -> primario UCSC;
        # mate del paciente -> segundo breakpoint UCSC.
        for ref in reference_index.get((local_chrom, mate_chrom), []):
            d_local = abs(local_pos - ref["primary_pos"])
            d_mate = abs(mate_pos - ref["mate_pos"])

            if d_local <= args.max_distance and d_mate <= args.max_distance:
                orientation_match = query_orientation == ref["orientation"]
                position_candidates.append({
                    "ref": ref,
                    "order": "direct",
                    "d_local": d_local,
                    "d_mate": d_mate,
                    "orientation_match": orientation_match,
                })

        # Swapped: mate del paciente -> primario UCSC;
        # local del paciente -> segundo breakpoint UCSC.
        for ref in reference_index.get((mate_chrom, local_chrom), []):
            d_local = abs(local_pos - ref["mate_pos"])
            d_mate = abs(mate_pos - ref["primary_pos"])

            if d_local <= args.max_distance and d_mate <= args.max_distance:
                orientation_match = (
                    query_orientation == tuple(reversed(ref["orientation"]))
                )
                position_candidates.append({
                    "ref": ref,
                    "order": "swapped",
                    "d_local": d_local,
                    "d_mate": d_mate,
                    "orientation_match": orientation_match,
                })

        # Evitar duplicar la misma combinación registro/orden.
        deduplicated = {}
        for candidate in position_candidates:
            ref = candidate["ref"]
            key = (
                ref["chrom_raw"], ref["start0"], ref["end0"],
                ref["name"], ref["second_bp"], candidate["order"],
            )
            deduplicated.setdefault(key, candidate)

        position_candidates = list(deduplicated.values())
        n_position = len(position_candidates)
        n_orientation = sum(
            c["orientation_match"] for c in position_candidates
        )

        if n_position == 0:
            status = "no_pair_within_distance"
        elif n_orientation == 0:
            status = "orientation_conflict"
        elif n_orientation == 1:
            status = "one_candidate"
        else:
            status = "multiple_candidates"

        summaries.append({
            **base,
            "STATUS": status,
            "N_POSITION_CANDIDATES": str(n_position),
            "N_ORIENTATION_MATCHES": str(n_orientation),
        })

        for candidate in position_candidates:
            ref = candidate["ref"]
            candidate_rows.append({
                "Samples_ID": sample,
                "SV_chrom": raw_chrom,
                "SV_start": raw_start,
                "SV_type": svtype,
                "ALT": alt,
                "SAMPLE_MATE_CHROM": mate_chrom,
                "SAMPLE_MATE_POS": str(mate_pos),
                "UCSC_ID": ref["name"],
                "UCSC_PRIMARY_COORD_BED": (
                    f'{ref["chrom_raw"]}:{ref["start0"]}-{ref["end0"]}'
                ),
                "UCSC_SECOND_BP": normalize_second_bp(ref["second_bp"]),
                "ENDPOINT_ORDER": candidate["order"],
                "DIST_LOCAL_BP": str(candidate["d_local"]),
                "DIST_MATE_BP": str(candidate["d_mate"]),
                "ORIENTATION_MATCH": str(candidate["orientation_match"]),
                "LR_AF": ref["AF"],
                "LR_AF_AFR": ref["AF_AFR"],
                "LR_AF_AMR": ref["AF_AMR"],
                "LR_AF_EAS": ref["AF_EAS"],
                "LR_AF_EUR": ref["AF_EUR"],
                "LR_AF_SAS": ref["AF_SAS"],
                "STATUS": status,
            })

    write_tsv(args.summary, summary_fields, summaries)
    write_tsv(args.candidates, candidate_fields, candidate_rows)

    print("Registros UCSC BND indexados:", n_reference)
    print("Llamadas BND/TRA únicas:", len(summaries))
    print("Parejas candidatas:", len(candidate_rows))
    print("Estados del resumen:")
    for status, count in sorted(Counter(r["STATUS"] for r in summaries).items()):
        print(f"  {status}: {count}")
    print("Resumen:", Path(args.summary).resolve())
    print("Candidatas:", Path(args.candidates).resolve())


if __name__ == "__main__":
    main()
