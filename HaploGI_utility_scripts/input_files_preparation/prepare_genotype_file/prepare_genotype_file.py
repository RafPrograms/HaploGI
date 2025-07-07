#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# prepare_genotype_file.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-27
#
# Usage:
#     python prepare_genotype_file.py --vcf <input.vcf or input.vcf.gz> \ 
#                                    [--output_dir <directory>]
#
# Required arguments:
#     --vcf            Path to the input VCF file (.vcf or .vcf.gz)
#
# Optional arguments:
#     --output_dir     Directory to save the output SNV genotypes file for HaploGI runs and log file (default: current directory)
#
# Description:
#     This script reads a VCF file and outputs a space-delimited SNV genotypes file.
#     Each variant is represented by a row, and each sample's alleles are coded as:
#         0 → 1
#         1 → 2
#         . (missing) → 0
#
#     The output is written to: <output_dir>/<input_basename>_SNV_genotypes.txt
#
# Example:
#     python prepare_genotype_file.py --vcf input_files/sample.vcf.gz \
#                                     --output_dir output_files_test/
#

import argparse
import gzip
import os
import logging


def setup_logging(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "prepare_genotype_file.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),  # Overwrite log file
            logging.StreamHandler()
        ]
    )


def parse_vcf(vcf_file):
    variants = []
    sample_ids = []

    open_func = gzip.open if vcf_file.endswith('.gz') else open

    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                sample_ids = line.strip().split('\t')[9:]
                continue

            fields = line.strip().split('\t')
            chrom = fields[0].replace("chr", "")
            pos = fields[1]
            variant_position = f"{chrom}:{pos}"
            sample_data = fields[9:]

            genotypes = []
            for sample in sample_data:
                genotype = sample.split(":")[0]
                alleles = genotype.replace('|', ' ').replace('/', ' ').split()

                converted = []
                for allele in alleles:
                    if allele == '1':
                        converted.append('2')
                    elif allele == '0':
                        converted.append('1')
                    elif allele == '.':
                        converted.append('0')
                genotypes.extend(converted)

            variants.append([variant_position] + genotypes)

    return variants, sample_ids


def write_output(variants, sample_ids, output_file):
    with open(output_file, 'w') as out:
        header = ["variant_position"]
        for sample in sample_ids:
            header.extend([sample, sample])
        out.write(" ".join(header) + "\n")

        for variant in variants:
            out.write(" ".join(variant) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Convert VCF to SNV genotype file used in HaploGI runns.")
    parser.add_argument('--vcf', required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument('--output_dir', default='.', help="Output directory (default: current directory)")
    args = parser.parse_args()

    setup_logging(args.output_dir)

    if not os.path.isfile(args.vcf):
        logging.error(f"VCF file not found: {args.vcf}")
        return

    base_name = os.path.splitext(os.path.basename(args.vcf))[0].replace('.vcf', '')
    output_file = os.path.join(args.output_dir, f"{base_name}_SNV_genotypes.txt")

    try:
        logging.info(f"Processing VCF: {args.vcf}")
        variants, sample_ids = parse_vcf(args.vcf)

        logging.info(f"Writing output to: {output_file}")
        write_output(variants, sample_ids, output_file)

        logging.info("Processing completed successfully.")
        print(f"Output written to: {output_file}")
    except Exception as e:
        logging.exception("An error occurred during processing.")


if __name__ == '__main__':
    main()
