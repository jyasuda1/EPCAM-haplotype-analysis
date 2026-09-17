#!/usr/bin/env python3
"""Convert a SHAPEIT .haps file to a headerless phased VCF."""

import argparse
import gzip
import sys
from pathlib import Path


def open_text(path: Path, mode: str):
    """Open plain text or gzip-compressed text based on the filename."""
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert SHAPEIT .haps records to headerless phased VCF records. "
            "The first five columns must be CHROM, ID, POS, REF, ALT; "
            "the remaining columns must be paired haplotype alleles."
        )
    )
    parser.add_argument("input", type=Path, help="Input .haps or .haps.gz file")
    parser.add_argument("output", type=Path, help="Output .vcf or .vcf.gz file")
    return parser.parse_args()


def normalize_allele(value: str, line_number: int) -> str:
    if value in {"0", "1"}:
        return value
    if value in {".", "?", "NA"}:
        return "."
    raise ValueError(
        f"line {line_number}: invalid haplotype allele {value!r}; "
        "expected 0, 1, ., ?, or NA"
    )


def convert(input_path: Path, output_path: Path) -> None:
    expected_haplotype_count = None

    with open_text(input_path, "rt") as source, open_text(output_path, "wt") as dest:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue

            fields = line.split()
            if len(fields) < 7:
                raise ValueError(
                    f"line {line_number}: at least 7 columns are required"
                )

            chrom, variant_id, pos, ref, alt = fields[:5]
            haplotypes = fields[5:]

            if len(haplotypes) % 2 != 0:
                raise ValueError(
                    f"line {line_number}: found an odd number "
                    f"({len(haplotypes)}) of haplotype columns"
                )

            if expected_haplotype_count is None:
                expected_haplotype_count = len(haplotypes)
            elif len(haplotypes) != expected_haplotype_count:
                raise ValueError(
                    f"line {line_number}: haplotype column count changed from "
                    f"{expected_haplotype_count} to {len(haplotypes)}"
                )

            genotypes = []
            for i in range(0, len(haplotypes), 2):
                allele1 = normalize_allele(haplotypes[i], line_number)
                allele2 = normalize_allele(haplotypes[i + 1], line_number)
                genotypes.append(f"{allele1}|{allele2}")

            vcf_fields = [
                chrom,
                pos,
                variant_id,
                ref,
                alt,
                ".",  # QUAL
                ".",  # FILTER
                ".",  # INFO
                "GT", # FORMAT
                *genotypes,
            ]
            dest.write("\t".join(vcf_fields) + "\n")


def main() -> int:
    args = parse_args()
    try:
        convert(args.input, args.output)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
