#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# create_phased_vcf.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-24

"""
Description:
    Augments a VCF (or VCF.GZ) file with phased genotypes from HaploGI haplotype data.
    Outputs a phased VCF file and a log file to the specified directory or current dir.

Usage:
    python create_phased_vcf.py \
        --haplotypes haplotype_sequences.txt \
        --positions SNV_genomic_positions.txt \
        --vcf input.vcf(.gz) \
        [--output_dir output_directory]
"""

__version__ = "1.0"

import argparse
import logging
import os
import gzip
from collections import defaultdict, OrderedDict

REMAP = {'1': '0', '2': '1', '7': '.', '3': '.', '0': '.'}

def remap_char(char):
    return REMAP.get(char, '.')

def setup_logging(log_path):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler()
        ]
    )

def read_sequence_file(filepath):
    region_data = []
    chrom = None
    raw_data = defaultdict(lambda: defaultdict(dict))

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            try:
                region, sample_id, sequence = line.split()
                chrom_part, pos_range = region.split(":")
                start, end = map(int, pos_range.split("-"))
                sample_prefix, hap_suffix = sample_id.split("_")
            except ValueError:
                logging.warning(f"Skipping malformed line in sequence file: {line}")
                continue

            if chrom is None:
                chrom = chrom_part
            raw_data[(start, end)][sample_prefix][hap_suffix] = sequence

    for (start, end), samples in raw_data.items():
        sample_data = {}
        for sample_id, haplotypes in samples.items():
            hap0 = haplotypes.get("0")
            hap1 = haplotypes.get("1")
            if not hap0 or not hap1 or len(hap0) != len(hap1):
                logging.warning(f"Incomplete or unequal-length haplotypes for {sample_id} in region {start}-{end}")
                continue
            combined_genotypes = [
                f"{remap_char(a)}|{remap_char(b)}"
                for a, b in zip(hap0, hap1)
            ]
            sample_data[sample_id] = combined_genotypes
        region_data.append((start, end, sample_data))

    return region_data

def read_position_blocks(filepath, block_size=20):
    position_blocks = OrderedDict()
    block = []
    block_id = 1

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            try:
                pos = int(line.split()[0])
            except ValueError:
                logging.warning(f"Skipping invalid line in position file: {line}")
                continue
            block.append(pos)
            if len(block) == block_size:
                position_blocks[block_id] = block
                block = []
                block_id += 1
        if block:
            position_blocks[block_id] = block

    return position_blocks

def build_combined_structure(position_blocks, region_data):
    combined = OrderedDict()
    for block_id, positions in position_blocks.items():
        start, end = positions[0], positions[-1]
        matched_samples = next(
            (samples for s, e, samples in region_data if s == start and e == end),
            None
        )
        if not matched_samples:
            logging.warning(f"No matching region for block {block_id} ({start}-{end})")
            continue

        for i, pos in enumerate(positions):
            if pos not in combined:
                combined[pos] = {}
            for sample_id, genotypes in matched_samples.items():
                combined[pos][sample_id] = genotypes[i]

    return combined

def open_vcf_file(vcf_path):
    if vcf_path.endswith('.gz'):
        return gzip.open(vcf_path, 'rt')
    return open(vcf_path, 'r')

def filter_and_augment_vcf(vcf_input_path, combined_structure, vcf_output_path):
    all_sample_ids = sorted(next(iter(combined_structure.values())).keys())
    vcf_lines = []
    header_written = False

    with open_vcf_file(vcf_input_path) as vcf_file:
        for line in vcf_file:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.strip().split('\t')
                original_sample_ids = header[9:]
                sample_to_index = {sid: idx + 9 for idx, sid in enumerate(original_sample_ids)}

                new_header = header[:9] + all_sample_ids
                vcf_lines.append("##fileformat=VCFv4.2")
                vcf_lines.append('\t'.join(new_header))
                header_written = True
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            try:
                pos = int(parts[1])
            except ValueError:
                continue

            if pos not in combined_structure:
                continue

            sample_genotypes = []
            for sid in all_sample_ids:
                value = combined_structure[pos].get(sid, ".|.")
                original = parts[sample_to_index[sid]] if sid in sample_to_index else None

                if value == ".|." and original:
                    if original in {".|.", "./."}:
                        sample_genotypes.append("./.")
                    elif original == "0/0":
                        sample_genotypes.append("0|0")
                    elif original == "1/1":
                        sample_genotypes.append("1|1")
                    else:
                        sample_genotypes.append(original)
                else:
                    sample_genotypes.append(value)

            new_line = parts[:9] + sample_genotypes
            vcf_lines.append('\t'.join(new_line))

    if not header_written:
        logging.warning("#CHROM header not found in input VCF. Inserting default header.")
        default_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + all_sample_ids
        vcf_lines.insert(0, '\t'.join(default_header))
        vcf_lines.insert(0, "##fileformat=VCFv4.2")

    with open(vcf_output_path, 'w') as out_file:
        out_file.write('\n'.join(vcf_lines) + '\n')

    logging.info(f"Final VCF written to: {vcf_output_path}")

def main():
    parser = argparse.ArgumentParser(description="Augment a VCF file with phased genotypes using HaploGI haplotypes.")
    parser.add_argument('--haplotypes', required=True, help='Path to HaploGI haplotype sequences file')
    parser.add_argument('--positions', required=True, help='Path to SNV genomic positions file')
    parser.add_argument('--vcf', required=True, help='Path to input VCF file (.vcf or .vcf.gz)')
    parser.add_argument('--output_dir', default='.', help='Directory to save output files (default: current directory)')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    output_basename = 'HaploGI_phased_WGS_data.vcf'
    output_path = os.path.join(args.output_dir, output_basename)
    log_path = os.path.join(args.output_dir, 'HaploGI_phased_WGS_data.log')

    setup_logging(log_path)
    logging.info("=== Creation of VCF with Phased WGS Data Started ===")

    region_data = read_sequence_file(args.haplotypes)
    position_blocks = read_position_blocks(args.positions)
    combined_structure = build_combined_structure(position_blocks, region_data)
    filter_and_augment_vcf(args.vcf, combined_structure, output_path)

    logging.info("=== Creation of VCF with Phased WGS Data Completed ===")

if __name__ == "__main__":
    main()
