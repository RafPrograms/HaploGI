#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# plot_haplotype_sharing.py
#
# Author: Rafael Nafikov
# Version: 1.0
# Date: 2025-06-27
#
# Description:
#     This Python utility plots the results of haplotype sharing evaluation
#     produced by HaploGI. It visualizes the extent of haplotype sharing
#     among subjects from the core set of cases, helping identify a shared
#     risk haplotype and its boundaries within a region of interest (ROI).
#
# Usage:
#     python plot_haplotype_sharing.py \
#         --pedigree <pedigree_file> \
#         --positions <positions_file> \
#         --haploshare <haploshare_file> \
#         --set <set_file> \
#         [--roi <roi_bounds>] \
#         [--x_start <start_window>] \
#         [--x_end <end_window>] \
#         [--outdir <output_directory>]
#
# Arguments:
#     --pedigree       Pedigree file (≥5 columns) used in HaploGI runs; must include subject ID (column 0) and phenotype (column 4).
#     --positions      SNV genomic positions file used in HaploGI runs.
#     --haploshare     Haplotype sharing matrix file: first column = genomic window numbers; remaining columns = subject IDs.
#     --set            Space-separated list of subject IDs forming the core set of cases.
#     --roi            (Optional) Region of interest (ROI) as cM positions, with the middle value indicating the linkage peak (e.g., "5.0-17.4-35.1").
#     --x_start        (Optional) Minimum genomic window number to include in the plot.
#     --x_end          (Optional) Maximum genomic window number to include in the plot.
#     --outdir         (Optional) Output directory for saving results (default: ./output).
#
# Output:
#     - genomic_windows.txt: Processed genomic windows summary
#     - gnuplot_data_page_*.dat: One or more data files for Gnuplot plotting
#     - gnuplot_output_page_*.pdf/.eps: Multi-page plots in PDF and EPS formats
#     - gnuplot_output_legend.pdf/.eps: Legend explaining subject classification
#
# Requirements:
#     - Python 3.x
#     - NumPy
#     - Gnuplot installed and accessible via system PATH
#
# License:
#     GNU General Public License v3.0 or later
#     See https://www.gnu.org/licenses/gpl-3.0.html
#

__version__ = "1.0"



import sys
import numpy as np
import io
import logging
import subprocess
import math
import os
import argparse
# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)


chrs_total = 0

def read_subject_phenotype(file_path):
    """Reads subject and phenotype columns from pedigree file, skipping comment lines."""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('**') or not line:
                continue  # Skip comment or empty lines
            parts = line.split()
            if len(parts) >= 5:
                try:
                    subject = int(parts[0])
                    phenotype = int(parts[4])
                    data.append([subject, phenotype])
                except ValueError:
                    continue  # Skip lines with invalid integer data
    return np.array(data)


def process_blocks(file_path):
    """Processes SNV genomic positions file into genomic window blocks of 20 lines and formats output."""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    a = int(parts[0])
                    b = float(parts[1])
                    data.append((a, b))
                except ValueError:
                    continue  # Skip malformed lines

    result = []
    num_blocks = len(data) // 20  # Only complete 20-line blocks
    for i in range(num_blocks):
        block = data[i*20:(i+1)*20]
        first_num = block[0][0]
        last_num = block[-1][0]
        id_col = f"{first_num}-{last_num}"
        seq_nums = [str(entry[0]) for entry in block]
        val_first = str(block[0][1])
        val_last = str(block[-1][1])
        block_line = [str(i + 1), id_col] + seq_nums + [val_first, val_last]
        result.append(block_line)
    return result


def write_genomic_windows(output_file, block_results):
    """Writes processed genomic window blocks to output file with header."""
    with open(output_file, 'w') as f:
        header = ["Genomic_window", "boundaries_bps"] + [f"{i}_bp" for i in range(1, 21)] + ["1_cM_pos", "20_cM_pos"]
        f.write('\t'.join(header) + '\n')
        for row in block_results:
            f.write('\t'.join(row) + '\n')


def make_unique_zero_based_except_first(names):
    """Keep the first name untouched. Deduplicate the rest with _0, _1, etc."""
    counts = {}
    unique_names = [names[0]]  # First column name unchanged
    for name in names[1:]:
        count = counts.get(name, 0)
        unique_name = f"{name}_{count}"
        unique_names.append(unique_name)
        counts[name] = count + 1
    return unique_names


def load_data_with_unique_names(filepath):
    with open(filepath, 'r') as f:
        header_line = f.readline()
        raw_headers = header_line.strip().split()
        unique_headers = make_unique_zero_based_except_first(raw_headers)
        rest_of_data = f.read()
    full_data = io.StringIO(" ".join(unique_headers) + "\n" + rest_of_data)
    dtype = [(name, 'i4') for name in unique_headers]
    data = np.genfromtxt(full_data, dtype=dtype, names=True)
    return data


def read_integer_list_from_file(file_path):
    """Read a single line of integers separated by spaces representing IDs for the core set of case."""
    with open(file_path, 'r') as f:
        line = f.readline().strip()
        if not line:
            return []
        return [int(x) for x in line.split()]


def build_subject_data_structure(subject_phenotype, integer_list):
    """
    Construct structured data:
    col1: subject IDs
    col2: phenotype
    col3: 1 if subject in integer_list else 0
    col4: based on phenotype and col3
    """
    data_list = []
    integer_set = set(integer_list)  # for fast lookup
    for entry in subject_phenotype:
        subj = entry[0]
        pheno = entry[1]
        present = 1 if subj in integer_set else 0
        if pheno == 0:
            col4 = 4
        elif pheno == 1:
            col4 = 3
        elif pheno == 2 and present == 0:
            col4 = 2
        elif pheno == 2 and present == 1:
            col4 = 1
        else:
            col4 = 0  # fallback if unexpected
        data_list.append([subj, pheno, present, col4])
    return np.array(data_list)


def save_gnuplot_chunk(filename, x_vals, series_list, col_names, offset=0):
    """
    Write a data file for gnuplot with header line containing column names:
    first col: x_vals name (usually genomic window)
    following cols: series names matching the columns (used for colors)
    """
    with open(filename, "w") as f:
        # Write header line
        header_line = [col_names[0]] + list(col_names[1:])
        f.write(" ".join(header_line) + "\n")

        for i in range(len(x_vals)):
            row = [str(x_vals[i])]
            for j, val in enumerate(series_list):
                v = val[i]
                # Write row index (1-based) if value is 1 or 2, else NaN
                row.append(str(offset + j + 1) if v in (1, 2) else "NaN")
            f.write(" ".join(row) + "\n")

def split_matrix_rows(matrix_data, windows_per_page=3200):
    """Split matrix data into multiple chunks, each with up to `windows_per_page` rows."""
    total_rows = len(matrix_data)
    chunks = []

    for start in range(0, total_rows, windows_per_page):
        end = start + windows_per_page
        chunk = matrix_data[start:end]
        chunks.append(chunk)
        
    return chunks


def split_and_save_gnuplot_pages_by_windows_and_columns(matrix_data, base_filename="gnuplot_data_page", max_windows=3200, max_subjects=16):
    """
    Splits matrix data into multiple files, each with up to `max_windows` rows
    and `max_subjects` columns (excluding x-axis column).
    Returns list of dicts with filename and plotting metadata.
    """
    names = matrix_data.dtype.names
    x_col = names[0]
    subject_names = names[1:]
    global chrs_total
    chrs_total = len (subject_names)
    total_windows = len(matrix_data)
    total_subjects = len(subject_names)

    outputs = []
    row_chunks = [
        (start, min(start + max_windows, total_windows))
        for start in range(0, total_windows, max_windows)
    ]
    col_chunks = [
        (start, min(start + max_subjects, total_subjects))
        for start in range(0, total_subjects, max_subjects)
    ]

    for row_idx, (row_start, row_end) in enumerate(row_chunks):
        row_slice = matrix_data[row_start:row_end]
        x_vals_chunk = row_slice[x_col]
        xmin = int(x_vals_chunk[0])
        xmax = int(x_vals_chunk[-1])

        for col_idx, (col_start, col_end) in enumerate(col_chunks):
            subject_chunk_names = subject_names[col_start:col_end]
            series_list = [row_slice[name] for name in subject_chunk_names]
            chunk_col_names = [x_col] + list(subject_chunk_names)

            file_idx = row_idx * len(col_chunks) + col_idx + 1
            filename = f"{base_filename}_{file_idx}.dat"
            save_gnuplot_chunk(filename, x_vals_chunk, series_list, chunk_col_names, offset=col_start)

            outputs.append({
                "filename": filename,
                "xmin": xmin,
                "xmax": xmax,
                "subject_offset": col_start,
                "subject_names": list(subject_chunk_names)
            })

    return outputs


def map_positions_to_window_numbers(positions, block_results):
    """
    Map each input position to the closest window number in block_results.
    block_results is a list of lists, with position at index 22 (23rd col),
    window number at index 0 (1st col).
    """
    mapped_windows = []
    for pos in positions:
        closest_diff = float('inf')
        closest_window = None
        for row in block_results:
            try:
                block_pos = float(row[22])
                window_num = int(row[0])
                diff = abs(block_pos - pos)
                if diff < closest_diff:
                    closest_diff = diff
                    closest_window = window_num
            except (ValueError, IndexError):
                continue
        mapped_windows.append(closest_window)
    return mapped_windows


def generate_gnuplot_script(data_file_info, build_subject_data, block_results, vertical_lines_windows=None, base_pdf_name="gnuplot_output"):
    color_map = {1: "rgb 'magenta'", 2: "rgb 'orange'", 3: "rgb 'blue'", 4: "rgb 'black'"}
    subject_color_lookup = {entry[0]: color_map.get(entry[3], "rgb 'grey'") for entry in build_subject_data}

    # Prepare x2tics (top) for genomic positions every 200 windows
    x2tics = []
    for row in block_results:
        try:
            window_num = int(row[0])
            genomic_pos = float(row[22])
            if window_num % 200 == 0:
                label = f"{genomic_pos:.1f}"
                x2tics.append(f"\"{label}\" {window_num}")
        except (ValueError, IndexError):
            continue

    tmp1 = 0  # Y-axis offset for pages
    script_lines = []

    for page_idx, info in enumerate(data_file_info):
        data_file = info['filename']
        xmin = info['xmin']
        xmax = info['xmax']
        subject_offset = info['subject_offset']
        col_names = info['subject_names']

        for terminal, ext in [("pdfcairo", "pdf"), ("epscairo", "eps")]:
            script_lines.append("reset")
            script_lines.append(f"set terminal {terminal} enhanced color font 'Arial,10'")
            script_lines.append("set termoption solid")  # <-- Added this line
            script_lines.append("set size 1,0.6")
            script_lines.append("set lmargin at screen 0.15")
            script_lines.append("set rmargin at screen 0.85")
            script_lines.append("set tmargin at screen 0.85")
            script_lines.append("set bmargin at screen 0.15")
            script_lines.append("set xlabel 'Genomic Window'")
            script_lines.append("set xtics 200 nomirror rotate by -45")
            script_lines.append("set xtics out")
            script_lines.append("set ytics out")
            script_lines.append("unset key")
            script_lines.append("set offset graph 0.0, graph 0.0, graph 0.0, graph 0.01")
            script_lines.append("set style data linespoints")
            script_lines.append("set label 1 'Genomic position, cM' at screen 0.5, screen 0.95 center font ',12'")
            script_lines.append("set x2range [:]")
            script_lines.append("set x2tics scale 0.5")
            if x2tics:
                script_lines.append("set x2tics rotate by 45 (" + ", ".join(x2tics) + ")")
            else:
                script_lines.append("set x2tics rotate by 45")
            script_lines.append("set x2tics out")

            if vertical_lines_windows:
                colors = ["rgb 'brown'", "rgb 'black'", "rgb 'brown'"]
                for i, win in enumerate(vertical_lines_windows):
                    if win is not None:
                        script_lines.append(f"set arrow {i+1} from {win}, graph 0 to {win}, graph 1 nohead lw 2 lc {colors[i]}")

            raw_labels = [label.rsplit('_', 1)[0] for label in col_names]
            y_labels = [f'"{label}" {idx+1+tmp1}' for idx, label in enumerate(raw_labels)]
            script_lines.append("set ytics (" + ", ".join(y_labels) + ")")
            start = 0.5 + tmp1
            end = len(y_labels) + tmp1 + 0.5
            script_lines.append(f"set yrange [{start}:{end}]")
            script_lines.append(f"set xrange [{xmin}:{xmax}]")

            plots = []
            for idx, col_name in enumerate(col_names, start=1):
                base_name = col_name.rsplit('_', 1)[0]
                try:
                    subject_id = int(base_name)
                except ValueError:
                    subject_id = None
                color = subject_color_lookup.get(subject_id, "rgb 'grey'")
                col_num = idx + 1
                # Changed lw from 25 to 7 here
                plots.append(f"'{data_file}' using 1:{col_num} with lines lt 1 lw 15 lc {color} title 'Series {subject_id}'")

            output_file = f"{base_pdf_name}_page_{page_idx + 1}.{ext}"
            script_lines.append(f"set output '{output_file}'")
            script_lines.append("plot " + ", \\\n     ".join(plots))
            script_lines.append("unset output\n")

        if tmp1 + len(col_names) < chrs_total:
            tmp1 += len(col_names)
        else:
            tmp1 = 0

    # Legend Page (PDF and EPS) unchanged, but also add set termoption solid here:

    legend_entries = {
        1: ("Case, part of the core group", "rgb 'magenta'"),
        2: ("Case, not part of the core group", "rgb 'orange'"),
        3: ("Control", "rgb 'blue'"),
        4: ("No phenotype data", "rgb 'black'")
    }

    for terminal, ext in [("pdfcairo", "pdf"), ("epscairo", "eps")]:
        script_lines.append("reset")
        script_lines.append(f"set terminal {terminal} enhanced color font 'Arial,10'")
        script_lines.append("set termoption solid")  # <-- Added here too
        script_lines.append(f"set output '{base_pdf_name}_legend.{ext}'")
        script_lines.append("unset border")
        script_lines.append("unset tics")
        script_lines.append("unset key")
        script_lines.append("set xrange [0:10]")
        script_lines.append(f"set yrange [0:{len(legend_entries) + 1}]")
        script_lines.append("set size 0.8,0.5")
        script_lines.append("set label 'Legend' at screen 0.5, screen 0.95 center font ',14'")
        script_lines.append("set multiplot")

        for i, (key, (label, color)) in enumerate(sorted(legend_entries.items()), start=1):
            y = len(legend_entries) - i + 1
            object_id = i
            label_id = 100 + i
            script_lines.append(f"set object {object_id} rect from 1,{y} to 2,{y + 0.5} fc {color} fs solid 1.0")
            script_lines.append(f"set label {label_id} '{label}' at 2.5,{y + 0.25} left")

        script_lines.append("plot NaN notitle")
        script_lines.append("unset multiplot")
        script_lines.append("unset output\n")

    return "\n".join(script_lines)


def run_gnuplot_script(script_text, script_filename="plot_script.gp"):
    with open(script_filename, 'w') as f:
        f.write(script_text)
    try:
        subprocess.run(["gnuplot", script_filename], check=True)
        logging.info(f"Gnuplot plotting done, see output PDFs.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Gnuplot failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate haplotype sharing plots using data files projueced by HaploGI."
    )
    parser.add_argument("--pedigree", required=True, help="Pedigree file with subject ID and phenotype info.")
    parser.add_argument("--positions", required=True, help="SNV genomic positions file.")
    parser.add_argument("--haploshare", required=True, help="Haplotype sharing data generated by HaploGI.")
    parser.add_argument("--set", required=True, help="Space-separated list of subject IDs from the core set of cases.")
    parser.add_argument("--roi", help="ROI as: left-mid-right cM positions with the middle value indicating the linkage peak (e.g., 5.0-17.4-35.1)")
    parser.add_argument("--x_start", type=int, help="Minimum genomic window number to plot.")
    parser.add_argument("--x_end", type=int, help="Maximum genomic window number to plot.")
    parser.add_argument("--outdir", default="output", help="Directory to store all output files (default: ./output)")

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.outdir, exist_ok=True)

    # Set up logging to file and console
    log_file = os.path.join(args.outdir, "plot_haplotype_sharing.log")
    file_handler = logging.FileHandler(log_file, mode='w')  # Overwrite each run
    console_handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logging.info(f"Logging initialized. Output directory: {args.outdir}")


    for file in [args.pedigree, args.positions, args.haploshare, args.set]:
        if not os.path.exists(file):
            logging.error(f"File not found: {file}")
            sys.exit(1)

    logging.info("Loading subject/phenotype data...")
    subject_data = read_subject_phenotype(args.pedigree)

    logging.info("Processing SNV genomic positions data...")
    block_results = process_blocks(args.positions)
    write_genomic_windows(os.path.join(args.outdir, "genomic_windows.txt"), block_results)

    logging.info("Loading subject IDs from the core set of cases list...")
    integer_list = read_integer_list_from_file(args.set)

    logging.info("Building subject data structure...")
    build_subject_data = build_subject_data_structure(subject_data, integer_list)

    logging.info("Loading matrix with haplotype sharing data...")
    matrix_data = load_data_with_unique_names(args.haploshare)

    if args.x_start is not None and args.x_end is not None:
        x_col = matrix_data.dtype.names[0]
        matrix_data = matrix_data[(matrix_data[x_col] >= args.x_start) & (matrix_data[x_col] <= args.x_end)]
        logging.info(f"Filtered matrix data to x-range: {args.x_start} - {args.x_end}")

    # Data output directory
    base_file_prefix = os.path.join(args.outdir, "gnuplot_data_page")
    data_file_info = split_and_save_gnuplot_pages_by_windows_and_columns(matrix_data, base_filename=base_file_prefix)

    vertical_windows = None
    if args.roi:
        try:
            roi_bounds = [float(x) for x in args.roi.split('-')]
            vertical_windows = map_positions_to_window_numbers(roi_bounds, block_results)
            if args.x_start is not None and args.x_end is not None and vertical_windows:
                vertical_windows = [
                    win if win is not None and args.x_start <= win <= args.x_end else None
                    for win in vertical_windows
                ]
            logging.info(f"Vertical positions mapped to window numbers: {vertical_windows}")
        except Exception as e:
            logging.error(f"Invalid vertical position format: {e}")

    logging.info("Generating gnuplot script...")
    base_pdf_name = os.path.join(args.outdir, "gnuplot_output")
    gnuplot_script_text = generate_gnuplot_script(
        data_file_info, build_subject_data, block_results, vertical_lines_windows=vertical_windows, base_pdf_name=base_pdf_name
    )

    script_path = os.path.join(args.outdir, "plot_script.gp")
    logging.info("Running gnuplot...")
    run_gnuplot_script(gnuplot_script_text, script_filename=script_path)


if __name__ == "__main__":
    main()
