# `create_genomic_windows.py`

![Python](https://img.shields.io/badge/python-3.x-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3.0-blue)

🔗 [`create_genomic_windows.py`](./create_genomic_windows.py)

**Version:** 1.0  
**Author:** Rafael Nafikov  
**Date:** 2025-06-20

---

## Table of Contents
- [Folder Structure](#folder-structure)
- [Overview](#overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Arguments](#arguments)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Example](#example)
- [License](#license)
- [Support](#support)

---

## Folder Structure

```
create_genomic_windows_file/
├── create_genomic_windows.py
├── input_files
│   └── chr16_SNV_genomic_positions_B38.txt
├── output_files
│   ├── create_genomic_windows.log
│   └── genomic_windows.txt
└── README.md

```

3 directories, 5 files

---

## Overview

`create_genomic_windows.py` is a Python utility designed to generate a genomic windows file for reference purposes.  
This file is particularly helpful when navigating the results of HaploGI runs.

The script takes as input an SNV genomic positions file — typically the same one prepared for HaploGI — and divides it into fixed-size blocks (default: 20 SNVs per block). Each block includes the genomic boundaries, the SNV positions within that block, and their corresponding genetic positions in centimorgans (cM).

---

## Requirements

- Python 3.x  
- No external dependencies (standard library only)

---

## Usage

Run the script from the command line with the required inputs:

```bash
python create_genomic_windows.py --positions <SNV_genomic_positions.txt> \
                                [--outdir <directory>]
```
 
---

## Arguments

| Argument      | Required | Description                                         |
|---------------|----------|-----------------------------------------------------|
| `--positions`   | Yes      | Path to the SNV genomic positions file used in HaploGI runs|
| `--outdir`   | No      | Directory where the output genomic windows file and log will be saved (default: ./output)|

---

## Input Files

### 🔷 `SNV_genomic_positions.txt`
```bash
1052701 3.767099
1052874 3.767696
1053095 3.768460
```

- Physical base pair position for each variant
- Genomic position of the variant in cM

---

## Output Files

### 🟢 `genomic_windows.txt`

```bash
Genomic_window	boundaries_bps	1_bp	2_bp	3_bp	4_bp	5_bp	6_bp	7_bp	8_bp	9_bp	10_bp	11_bp	12_bp	13_bp	14_bp	15_bp	16_bp	17_bp	18_bp	19_bp	20_bp	1_cM_pos	20_cM_pos
1	1052701-1055604	1052701	1052874	1053095	1053154	1053316	1053758	1053890	1054244	1054249	1054409	1054446	1054491	1054508	1054606	1054612	1054759	1055201	1055294	1055431	1055604	3.767099	3.77713
2	1055638-1059301	1055638	1055760	1055783	1055819	1055825	1055873	1055899	1055919	1055937	1055966	1056064	1057345	1057865	1058373	1058485	1058669	1058884	1059044	1059178	1059301	3.777248	3.789907
```
Where each row represents a genomic window, and the columns are as follows:

- Genomic_window: Sequential genomic window number.

- boundaries_bps: The genomic positions marking the start and end of the window (1st variant base pair position to the 20th variant base pair position).

- 1_bp to 20_bp: Base pair positions for 20 variants comprising the genomic window.

- 1_cM_pos and 20_cM_pos: The genomic position in centimorgan (cM) for the 1st and 20th variants within the genomic window.

### 🟢 `create_genomic_windows.log`
A log summary of the processing steps, including any warnings and errors encountered during execution.

---

## Example

Unzip the `input_files` directory before running the script. Use the included example inputs to run:

```bash
python create_genomic_windows.py --positions ./input_files/chr16_SNV_genomic_positions_B38.txt \
                                 --outdir ./output_files_test
```

📌 Note: If you specify `--output_dir`, make sure that directory exists before running the script:
```bash
 mkdir -p ./output_files_test
```

---

## License

This project is licensed under the **GNU General Public License v3.0 or later** (`GPL-3.0-or-later`).  
It uses the SPDX license identifier `GPL-3.0-or-later` in source files for standardized license declaration.

You can view the full license text at:  
https://www.gnu.org/licenses/gpl-3.0.html

Unless required by applicable law or agreed to in writing, the software is distributed  
**"AS IS"**, without warranties or conditions of any kind, either express or implied.

See the [LICENSE](https://github.com/RafPrograms/HaploGI/blob/main/LICENSE) file in this repository for complete details.

---

## Support

For questions, bug reports, or suggestions, please contact:  
📧 your.email@example.com  
🔗 [GitHub Issues](https://github.com/yourusername/yourrepo/issues)

