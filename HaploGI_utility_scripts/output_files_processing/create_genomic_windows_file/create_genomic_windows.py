#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# create_genomic_windows.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-28
#
# Usage:
#     python create_genomic_windows.py --positions <positions.txt> \
#                                     [--outdir <directory>]
#
# Required arguments:
#     --positions       Input file with SNV genomic positions (tab-delimited, with at least 2 columns: positions in base pairs and cM)
#
# Optional arguments:
#     --outdir          Output directory for results and logs (default: ./output)
#
# Description:
#     This script processes a file of SNV genomic positions and divides them into fixed-size blocks
#     (default block size: 20). It outputs a genomic windows file where each row corresponds to a block,
#     listing block boundaries, SNV positions, and cM positions of the first and last SNVs.
#
#     Each genomic window includes:
#         - Genomic_window number
#         - Block boundaries in base pairs
#         - 20 SNV base pair positions
#         - Start and end cM positions
#
# Output:
#     - <outdir>/genomic_windows.txt  — Tab-delimited file with genomic window definitions
#     - <outdir>/create_genomic_windows.log  — Processing log
#
# Example:
#     python create_genomic_windows.py --positions data/positions_chr16.txt \
#                                      --outdir results/


__version__ = "1.0"

import sys
import os
import argparse
import logging

# Constants
BLOCK_SIZE = 20

# Logger setup
logger = logging.getLogger("GenomicPlotGenerator")
logger.setLevel(logging.INFO)

def setup_logging(log_file):
    """Sets up logging to both file and console."""
    logger.handlers.clear()
    file_handler = logging.FileHandler(log_file, mode='w')
    console_handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

def process_blocks(file_path):
    """Processes genomic position file into fixed-size blocks."""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    pos = int(parts[0])
                    cm = float(parts[1])
                    data.append((pos, cm))
                except ValueError:
                    continue  # Skip lines that can't be parsed

    result = []
    num_blocks = len(data) // BLOCK_SIZE
    for i in range(num_blocks):
        block = data[i * BLOCK_SIZE:(i + 1) * BLOCK_SIZE]
        id_col = f"{block[0][0]}-{block[-1][0]}"
        block_line = [str(i + 1), id_col] + [str(entry[0]) for entry in block] + [str(block[0][1]), str(block[-1][1])]
        result.append(block_line)

    logger.info(f"Processed {num_blocks} full blocks (of {BLOCK_SIZE} lines each).")
    return result

def write_genomic_windows(output_file, blocks):
    """Writes the processed genomic windows to a file."""
    header = (
        ["Genomic_window", "boundaries_bps"] +
        [f"{i}_bp" for i in range(1, BLOCK_SIZE + 1)] +
        ["1_cM_pos", "20_cM_pos"]
    )
    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in blocks:
            f.write('\t'.join(row) + '\n')

    logger.info(f"Output written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate genomic windows file from SNV positions.")
    parser.add_argument("--positions", required=True, help="Input file with SNV genomic positions.")
    parser.add_argument("--outdir", default="output", help="Output directory (default: ./output)")
    args = parser.parse_args()

    # Prepare output directory
    os.makedirs(args.outdir, exist_ok=True)
    log_path = os.path.join(args.outdir, "create_genomic_windows.log")
    setup_logging(log_path)

    if not os.path.isfile(args.positions):
        logger.error(f"Input file not found: {args.positions}")
        sys.exit(1)

    logger.info("Starting genomic block processing...")
    blocks = process_blocks(args.positions)

    output_file = os.path.join(args.outdir, "genomic_windows.txt")
    write_genomic_windows(output_file, blocks)

    logger.info("Processing complete.")

if __name__ == "__main__":
    main()
