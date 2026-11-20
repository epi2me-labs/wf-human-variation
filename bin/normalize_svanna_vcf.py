#!/usr/bin/env python3
import gzip
import re
import sys

src, dst = sys.argv[1:3]

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

contigs = []
with open_vcf(src) as handle:
    for line in handle:
        if not line.startswith("#"):
            break
        match = re.match(r"##contig=<ID=([^,>]+)", line)
        if match:
            contigs.append(match.group(1))

by_lower = {}
for contig in contigs:
    by_lower.setdefault(contig.casefold(), []).append(contig)

pattern = re.compile(r"([\[\]])([^:\[\]]+):([0-9]+)\1")
stats = {"records": 0, "mates": 0}
unknown = set()

with open_vcf(src) as fin, open(dst, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 8:
            fout.write(line)
            continue

        original_alt = fields[4]

        def normalize(match):
            mate = match.group(2)
            candidates = by_lower.get(mate.casefold(), [])

            if len(candidates) == 1:
                canonical = candidates[0]
                if canonical != mate:
                    stats["mates"] += 1
                return (
                    match.group(1) + canonical + ":"
                    + match.group(3) + match.group(1)
                )

            unknown.add(mate)
            return match.group(0)

        fields[4] = ",".join(
            pattern.sub(normalize, alt)
            for alt in fields[4].split(",")
        )

        if fields[4] != original_alt:
            stats["records"] += 1

        fout.write("\t".join(fields) + "\n")

print(f"Variantes BND modificadas: {stats['records']}")
print(f"Nombres de contig normalizados: {stats['mates']}")
print(f"VCF creado: {dst}")

if unknown:
    print("Contigs remotos sin coincidencia única:")
    for contig in sorted(unknown):
        print("  " + contig)
