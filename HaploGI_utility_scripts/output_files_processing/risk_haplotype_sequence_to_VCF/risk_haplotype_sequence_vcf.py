#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later

# risk_haplotype_sequence_vcf.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-27
#
# Description:
#     Processes a HaploGI-generated risk (shared) haplotype file, SNV genomic positions file,
#     and VCF used to generate a genotype file for HaploGI runs to produce:
#       - A genomic windows file (for your reference)
#       - A minimal VCF with phased risk haplotype sequences
#
# Usage:
#     python risk_haplotype_sequence_vcf.py \
#         --start_pos <start_position> \
#         --end_pos <end_position> \
#         --haplotypes <risk_haplotype_file> \
#         --vcf <input_vcf_file(.vcf or .vcf.gz)> \
#         --positions <position_blocks_file> \
#         [--output_dir <folder_path>]

__version__ = "1.0"

import argparse
import logging
import os
import gzip


def parse_and_filter_risk_haplotype_file(filepath, user_start, user_end):
    parsed_data = []
    found_start = False

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            region, value = parts[0], parts[1]

            try:
                coords = region.split(':')[1]
                start, end = map(int, coords.split('-'))

                if not found_start:
                    if start == user_start:
                        found_start = True
                        parsed_data.append((start, end, value))
                        if end == user_end:
                            break
                else:
                    parsed_data.append((start, end, value))
                    if end == user_end:
                        break

            except (IndexError, ValueError):
                logging.warning(f"Skipping malformed line: {line}")
                continue

    return parsed_data


def read_positions_file_and_create_blocks(pos_path, block_size=20):
    positions = []
    with open(pos_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    pos1 = int(parts[0])
                    pos2 = float(parts[1])
                    positions.append((pos1, pos2))
                except ValueError:
                    logging.warning(f"Skipping malformed line: {line}")

    blocks = []
    for i in range(0, len(positions) - block_size + 1, block_size):
        block = positions[i:i + block_size]
        if len(block) == block_size:
            row = [i // block_size + 1,
                   f"{block[0][0]}-{block[-1][0]}"] + \
                  [p[0] for p in block] + \
                  [block[0][1], block[-1][1]]
            blocks.append(row)

    return blocks


def write_blocks_to_file(blocks, output_path):
    with open(output_path, "w") as out:
        header = ["genomic_window", "boundaries"] + \
                 [f"var{i+1}" for i in range(20)] + \
                 ["var1_cM_pos", "var20_cM_pos"]
        out.write("\t".join(header) + "\n")
        for block in blocks:
            out.write("\t".join(map(str, block)) + "\n")
    logging.info(f"Genomic windows data written to: {output_path}")



def read_vcf_file(vcf_path):
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt" if vcf_path.endswith(".gz") else "r"

    with open_fn(vcf_path, mode) as vcf:
        return vcf.readlines()


def filter_blocks_by_haplotype_positions(blocks, haplotypes):
    filtered_blocks = []

    for idx, (haplo_start, haplo_end, _) in enumerate(haplotypes, start=1):
        match = next(
            (b for b in blocks if int(b[2]) == haplo_start and int(b[21]) == haplo_end),
            None
        )
        if match:
            filtered_blocks.append(match)
        else:
            logging.warning(f"No matching block for haplotype {idx} ({haplo_start}-{haplo_end})")

    if len(filtered_blocks) != len(haplotypes):
        logging.warning("Filtered blocks count does not match haplotypes count.")

    return filtered_blocks


def extract_window_and_haplotype_matrix(blocks, haplotypes):
    if len(blocks) != len(haplotypes):
        logging.warning("Mismatch between block and haplotype count.")
        return []

    allele_map = {'1': '0|.', '2': '1|.', '0': './.'}
    matrix = []

    for idx, (block, (_, _, haplo_str)) in enumerate(zip(blocks, haplotypes), 1):
        positions = block[2:22]
        alleles = list(haplo_str.strip())

        if len(alleles) != 20:
            logging.warning(f"Skipping HAP{idx}: invalid length")
            continue

        mapped = [allele_map.get(a, a) for a in alleles]
        matrix.append((f"HAP{idx}", list(zip(positions, mapped))))

    return matrix


def write_minimal_vcf(vcf_lines, matrix, output_vcf):
    genotype_by_pos = {int(pos): allele for _, row in matrix for pos, allele in row}

    with open(output_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##phasing=partial\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trisk_haplotype\n")

        for line in vcf_lines:
            if line.startswith("#"):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                logging.warning(f"Skipping malformed VCF line: {line.strip()}")
                continue

            chrom, pos, vid, ref, alt = fields[:5]
            pos = int(pos)

            if pos in genotype_by_pos:
                gt = genotype_by_pos[pos]
                out.write(f"{chrom}\t{pos}\t{vid or '.'}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

    logging.info(f"Risk haplotype VCF written to: {output_vcf}")


def main():
    parser = argparse.ArgumentParser(description="Process haplotypes and genomic data.")
    parser.add_argument("--start_pos", type=int, required=True, help="Start position (bp)")
    parser.add_argument("--end_pos", type=int, required=True, help="End position (bp)")
    parser.add_argument("--haplotypes", required=True, help="Risk haplotype file path")
    parser.add_argument("--vcf", required=True, help="VCF file path (.vcf or .vcf.gz)")
    parser.add_argument("--positions", required=True, help="Positions file path")
    parser.add_argument("--output_dir", default=".", help="Directory to save output files")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Set up logging to file and console
    log_file = os.path.join(args.output_dir, "risk_haplotype_sequence_vcf.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

    # Validate input files
    for path in [args.haplotypes, args.vcf, args.positions]:
        if not os.path.isfile(path):
            logging.error(f"Input file not found: {path}")
            raise FileNotFoundError(f"Input file not found: {path}")

    haplotypes = parse_and_filter_risk_haplotype_file(args.haplotypes, args.start_pos, args.end_pos)
    logging.info(f"Loaded {len(haplotypes)} lines from shared haplotype file.")

    vcf_lines = read_vcf_file(args.vcf)
    logging.info(f"VCF file loaded: {len(vcf_lines)} lines.")

    blocks = read_positions_file_and_create_blocks(args.positions)
    logging.info(f"Generated {len(blocks)} genomic windows.")

    windows_file = os.path.join(args.output_dir, "genomic_windows.txt")
    vcf_file = os.path.join(args.output_dir, "risk_haplotype.vcf")

    write_blocks_to_file(blocks, output_path=windows_file)

    filtered_blocks = filter_blocks_by_haplotype_positions(blocks, haplotypes)
    logging.info(f"{len(filtered_blocks)} genomic windows within specified boundaries had haplotype data.")

    matrix = extract_window_and_haplotype_matrix(filtered_blocks, haplotypes)

    if matrix:
        write_minimal_vcf(vcf_lines, matrix, output_vcf=vcf_file)


if __name__ == "__main__":
    main()
