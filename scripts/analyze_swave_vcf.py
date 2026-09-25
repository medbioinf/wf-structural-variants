#!/usr/bin/env python3
"""
Analyzes a Swave-annotated VCF (with REPEAT_RATIO from TRF and gene
annotations from ANNOVAR) and reports SSV/CSV counts, singleton and
allele-frequency tier breakdowns, and repeat/exon-disruption filtering,
mirroring the reporting style used in the original Swave paper.
"""

import argparse
import gzip

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str):
    """Parses a VCF INFO string into a dict, keeping flag-only keys as True."""
    info = {}
    for field in info_str.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            info[key] = value
        else:
            info[field] = True
    return info


def is_exon_disrupting(info):
    """
    Mirrors the paper's approach: "SVs intersecting 'exons' or 'splicing'
    annotations were flagged for downstream analysis."
    """
    func = info.get("Func.refGene", "")
    return any(tag in func for tag in ("exonic", "splicing"))


def classify_af(af):
    if af < 0.01:
        return "rare"
    elif af < 0.05:
        return "infrequent"
    else:
        return "frequent"


def analyze_vcf(vcf_path, repeat_threshold=0.8):
    stats = {
        "total": 0,
        "ssv": 0,
        "csv": 0,
        "ssv_singleton": 0,
        "csv_singleton": 0,
        "ssv_af": {"rare": 0, "infrequent": 0, "frequent": 0},
        "csv_af": {"rare": 0, "infrequent": 0, "frequent": 0},
        "ssv_highly_repetitive": 0,
        "csv_highly_repetitive": 0,
        "ssv_exon_disrupting": 0,
        "csv_exon_disrupting": 0,
        "ssv_singleton_exon_disrupting": 0,
        "csv_singleton_exon_disrupting": 0,
    }

    with open_vcf(vcf_path) as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            info_str = parts[7]
            info = parse_info(info_str)

            svtype = info.get("SVTYPE", "")
            is_csv = "+" in svtype
            category = "csv" if is_csv else "ssv"

            stats["total"] += 1
            stats[category] += 1

            ac = int(info.get("AC", 0))
            is_singleton = ac == 1
            if is_singleton:
                stats[f"{category}_singleton"] += 1

            af = float(info.get("AF", 0.0))
            af_tier = classify_af(af)
            stats[f"{category}_af"][af_tier] += 1

            repeat_ratio = info.get("REPEAT_RATIO")
            if repeat_ratio is not None and float(repeat_ratio) > repeat_threshold:
                stats[f"{category}_highly_repetitive"] += 1

            if is_exon_disrupting(info):
                stats[f"{category}_exon_disrupting"] += 1
                if is_singleton:
                    stats[f"{category}_singleton_exon_disrupting"] += 1

    return stats


def pct(part, whole):
    if whole == 0:
        return 0.0
    return round(100 * part / whole, 1)


def print_report(stats, repeat_threshold):
    total_ssv = stats["ssv"]
    total_csv = stats["csv"]
    total = stats["total"]

    print(f"\nTotal variants: {total}")
    print(f"  SSVs: {total_ssv} ({pct(total_ssv, total)}%)")
    print(f"  CSVs: {total_csv} ({pct(total_csv, total)}%)")

    print("\n--- Allele frequency tiers ---")
    for category, label in [("ssv_af", "SSV"), ("csv_af", "CSV")]:
        af = stats[category]
        cat_total = stats[category.split("_")[0]]
        print(f"\n{label} (n={cat_total}):")
        print(f"  Rare       (AF < 1%): {af['rare']} ({pct(af['rare'], cat_total)}%)")
        print(f"  Infrequent (1-5%):    {af['infrequent']} ({pct(af['infrequent'], cat_total)}%)")
        print(f"  Frequent   (> 5%):    {af['frequent']} ({pct(af['frequent'], cat_total)}%)")

    print("\n--- Singletons (AC=1) ---")
    print(f"  SSV singletons: {stats['ssv_singleton']} ({pct(stats['ssv_singleton'], total_ssv)}% of SSVs)")
    print(f"  CSV singletons: {stats['csv_singleton']} ({pct(stats['csv_singleton'], total_csv)}% of CSVs)")

    print(f"\n--- Repeat annotation (threshold: REPEAT_RATIO > {repeat_threshold}) ---")
    ssv_rep = stats["ssv_highly_repetitive"]
    csv_rep = stats["csv_highly_repetitive"]
    print(f"  SSVs flagged as highly repetitive: {ssv_rep} ({pct(ssv_rep, total_ssv)}% of SSVs)")
    print(f"  CSVs flagged as highly repetitive: {csv_rep} ({pct(csv_rep, total_csv)}% of CSVs)")
    print(f"  SSVs remaining after repeat filter: {total_ssv - ssv_rep}")
    print(f"  CSVs remaining after repeat filter: {total_csv - csv_rep}")

    print("\n--- Exon/splicing disruption (ANNOVAR Func.refGene) ---")
    print(f"  SSVs exon-disrupting: {stats['ssv_exon_disrupting']} ({pct(stats['ssv_exon_disrupting'], total_ssv)}% of SSVs)")
    print(f"  CSVs exon-disrupting: {stats['csv_exon_disrupting']} ({pct(stats['csv_exon_disrupting'], total_csv)}% of CSVs)")

    print("\n--- Singleton AND exon-disrupting (candidate pathogenic) ---")
    print(f"  SSV: {stats['ssv_singleton_exon_disrupting']}")
    print(f"  CSV: {stats['csv_singleton_exon_disrupting']}")


def main():
    parser = argparse.ArgumentParser(description="Analyze a Swave VCF for SSV/CSV, AF, singleton, repeat, and exon-disruption statistics.")
    parser.add_argument("--vcf", required=True, help="Path to the (optionally gzipped) Swave VCF to analyze.")
    parser.add_argument("--repeat_threshold", type=float, default=0.8, help="REPEAT_RATIO threshold above which a variant is considered highly repetitive (default: 0.8, matching the original paper).")

    args = parser.parse_args()

    stats = analyze_vcf(args.vcf, repeat_threshold=args.repeat_threshold)
    print_report(stats, args.repeat_threshold)


if __name__ == "__main__":
    main()