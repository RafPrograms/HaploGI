#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# decrease_number_of_MI_iterations.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-27
#
# Description:
#    Decrease the number of iterations of Meiosis Indicators (MI) in the MI input file 
#    for use in HaploGI runs. During runtime, the user will be prompted to specify 
#    how many MI iterations to retain in the output MI file.
#
# Usage:
#     python decrease_number_of_MI_iterations.py --MI matrix.txt \
#                                               [--output_dir results]
#
# Arguments:
#     --MI           Input haplotype sharing matrix file (required)
#     --output_dir   Optional output directory (default: current directory)
#
# Output:
#     - selected_<N>_blocks_every_<K>_th_block.txt
#     - extract_mi_blocks.log
#
# Requirements:
#     - Python 3.x
#
# License:
#     GNU General Public License v3.0 or later
#     See https://www.gnu.org/licenses/gpl-3.0.html
#

__version__ = "1.0"

import argparse
import logging
import os

def setup_logging(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "decrease_number_of_MI_iterations.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def extract_blocks(lines, col1_val):
    blocks = []
    i = 0
    while i < len(lines) - 1:
        line1 = lines[i].strip().split()
        line2 = lines[i + 1].strip().split()
        if line1 and line2 and line1[0] == col1_val and line2[0] == col1_val:
            block_start = i
            i += 2
            block = lines[block_start:block_start+2]
            while i < len(lines) - 1:
                next1 = lines[i].strip().split()
                next2 = lines[i + 1].strip().split()
                if next1 and next2 and next1[0] == col1_val and next2[0] == col1_val:
                    break
                block.append(lines[i])
                i += 1
            blocks.append(block)
        else:
            i += 1
    return blocks

def main():
    parser = argparse.ArgumentParser(description="Extract every N-th block from a haplotype matrix.")
    parser.add_argument("--MI", required=True, help="Path to Meiosis Indicators (MI) file")
    parser.add_argument("--output_dir", default=".", help="Directory to save output and log (default: current directory)")
    args = parser.parse_args()

    setup_logging(args.output_dir)

    if not os.path.isfile(args.MI):
        logging.error(f"File not found: {args.MI}")
        return

    try:
        with open(args.MI, 'r') as f:
            lines = f.readlines()

        if not lines:
            logging.warning("The input file is empty.")
            return

        first_line = lines[0].strip().split()
        if not first_line:
            logging.warning("First line is empty or malformed.")
            return
        col1_val = first_line[0]

        blocks = extract_blocks(lines, col1_val)
        total_blocks = len(blocks)

        if total_blocks == 0:
            logging.warning("No valid blocks found in the input file.")
            return

        logging.info(f"Total blocks found: {total_blocks}")
        print(f"Total blocks found: {total_blocks}")

        while True:
            try:
                desired_count = int(input("Enter number of blocks to keep: "))
                if desired_count <= 0:
                    print("Please enter a positive number.")
                    continue
                if desired_count > total_blocks:
                    print(f"Only {total_blocks} blocks available. Please enter a number ≤ {total_blocks}.")
                    continue

                step = total_blocks // desired_count
                selected_blocks = []
                i = 0
                while len(selected_blocks) < desired_count and i < total_blocks:
                    selected_blocks.append(blocks[i])
                    i += step

                output_file = os.path.join(args.output_dir, f"selected_{desired_count}_blocks_every_{step}_th_block.mi")
                with open(output_file, 'w') as out:
                    for block in selected_blocks:
                        out.writelines(block)

                logging.info(f"Selected {len(selected_blocks)} blocks (every {step}-th block).")
                logging.info(f"Output saved to: {output_file}")
                print(f"Selected {len(selected_blocks)} blocks (every {step}-th block).")
                print(f"Output saved to: {output_file}")
                break

            except ValueError:
                print("Please enter a valid integer.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == '__main__':
    main()
