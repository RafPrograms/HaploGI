#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# risk_alleles_variants.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-20
#
# Description:
#     This script identifies genetic variants where alleles are uniquely present
#     on a disease-associated (risk) haplotype. Specifically, it detects variants
#     that are consistently heterozygous (1/2 or 2/1) in all individuals from a core
#     set of cases and homozygous (1/1 or 2/2) in all other included individuals.
#     Controls known to carry the risk haplotype can be excluded to refine specificity.
#     A VCF file (plain or .gz) is used to enrich the output with variant metadata.
#
# Usage:
#     python risk_alleles_variants.py --start_pos <int> \
#                                     --end_pos <int> \
#                                     --genotypes <genotype_file> \
#                                     --set <core_cases_file> \
#                                     --vcf <vcf_file(.vcf or .vcf.gz)> \
#                                     [--controls <shared_haplotype_controls_file>] \
#                                     [--output_dir <directory>]
#
# Output:
#     A tab-delimited text file named `risk_alleles_variants.txt` and a log file
#     `risk_alleles_variants.log`, both written to the output directory.
#
# Requirements:
#     - Python 3.6+
#     - pandas
#
# License:
#     GNU General Public License v3.0 or later
#     https://www.gnu.org/licenses/gpl-3.0.html
#

__version__ = "1.0"

import pandas as pd
import os
import argparse
import logging
import gzip

def setup_logging(output_dir):
    log_path = os.path.join(output_dir, "risk_alleles_variants.log")
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )

def read_list_file(file_path):
    if not file_path or not os.path.isfile(file_path):
        logging.warning(f"List file not found: {file_path}")
        return []
    with open(file_path, 'r') as f:
        return f.readline().strip().split()

def extract_position_column(df, col_name):
    return df[col_name].astype(str).str.split(':', n=1).str[-1].astype(int)

def get_genotype_column_pairs(df, sample_names):
    genotype_pairs = {}
    for name in sample_names:
        idx1 = df.columns.get_loc(name)
        idx2 = idx1 + 1
        genotype_pairs[name] = list(zip(
            df.iloc[:, idx1].astype(int),
            df.iloc[:, idx2].astype(int)
        ))
    return genotype_pairs

def get_samples_to_include(df, excluded_samples):
    return [
        df.columns[i]
        for i in range(1, df.shape[1], 2)
        if df.columns[i] not in excluded_samples
    ]

def all_heterozygous(pairs):
    return all(pair in {(1, 2), (2, 1)} for pair in pairs)

def all_homozygous(pairs):
    return all(pair in {(1, 1), (2, 2)} for pair in pairs)

def infer_risk_allele(pairs):
    if all(pair == (2, 2) for pair in pairs):
        return 0
    elif all(pair == (1, 1) for pair in pairs):
        return 1
    return 0

def identify_risk_variants(start_pos, end_pos, genotypes_file, case_file, control_file=None):
    logging.info("Loading genotype and sample data...")

    case_samples = set(read_list_file(case_file))
    control_samples = set(read_list_file(control_file)) if control_file else set()

    try:
        df = pd.read_csv(genotypes_file, sep=r"\s+")
    except Exception as e:
        logging.error(f"Could not read genotype file: {e}")
        return []

    position_col = df.columns[0]
    positions = extract_position_column(df, position_col)

    region_mask = (positions >= start_pos) & (positions <= end_pos)
    df = df[region_mask].reset_index(drop=True)
    positions = positions[region_mask].reset_index(drop=True)

    included_samples = get_samples_to_include(df, excluded_samples=control_samples)
    genotype_data = get_genotype_column_pairs(df, included_samples)

    results = []
    for idx in range(len(df)):
        case_genotypes = [genotype_data[col][idx] for col in included_samples if col in case_samples]
        other_genotypes = [genotype_data[col][idx] for col in included_samples if col not in case_samples]

        if not case_genotypes or not other_genotypes:
            continue
        if not all_heterozygous(case_genotypes):
            continue
        if not all_homozygous(other_genotypes):
            continue

        inferred_allele = infer_risk_allele(other_genotypes)
        results.append((str(positions[idx]), inferred_allele))

    logging.info(f"Found {len(results)} candidate variants with risk-specific alleles.")
    return results

def parse_vcf(vcf_file):
    vcf_info = {}
    try:
        open_fn = gzip.open if vcf_file.endswith('.gz') else open
        with open_fn(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    pos = str(int(parts[1]))
                    vcf_info[pos] = (parts[2], parts[3], parts[4])
    except Exception as e:
        logging.error(f"Error parsing VCF file: {e}")
    return vcf_info

def write_output(results, vcf_data, output_path):
    try:
        with open(output_path, 'w') as f:
            f.write("Position\tRisk_Allele(0=REF;1=ALT)\tID\tREF\tALT\n")
            count = 0
            for pos, allele in results:
                if pos in vcf_data:
                    id_, ref, alt = vcf_data[pos]
                    f.write(f"{pos}\t{allele}\t{id_}\t{ref}\t{alt}\n")
                    count += 1

            if count == 0:
                f.write("# No matching variants found in VCF.\n")
                logging.warning("No overlapping positions found between results and VCF.")

        logging.info(f"Results written to: {output_path}")
    except Exception as e:
        logging.error(f"Failed to write results: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Identify variants with alleles unique to a disease-associated haplotype."
    )
    parser.add_argument("--start_pos", type=int, required=True, help="Start position (inclusive)")
    parser.add_argument("--end_pos", type=int, required=True, help="End position (inclusive)")
    parser.add_argument("--genotypes", required=True, help="Genotype input file (space-separated format)")
    parser.add_argument("--set", required=True, help="File listing case sample names (core set)")
    parser.add_argument("--vcf", required=True, help="VCF file to provide variant metadata (.vcf or .vcf.gz)")
    parser.add_argument("--controls", help="Optional file listing controls sharing the risk haplotype")
    parser.add_argument("--output_dir", default=".", help="Directory to write output files")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    logging.info("Starting risk allele variant detection...")

    results = identify_risk_variants(
        args.start_pos,
        args.end_pos,
        args.genotypes,
        args.set,
        args.controls
    )

    vcf_data = parse_vcf(args.vcf)

    output_file = os.path.join(args.output_dir, "risk_alleles_variants.txt")
    write_output(results, vcf_data, output_file)

    logging.info("Pipeline completed successfully.")

if __name__ == "__main__":
    main()
