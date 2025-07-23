# `prepare_genotype_file.py`

![Python](https://img.shields.io/badge/python-3.x-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3.0-blue)

🔗 [`prepare_genotype_file.py`](./prepare_genotype_file.py)

**Version:** 1.0  
**Author:** Rafael Nafikov  
**Date:** 2025-06-27

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
prepare_genotype_file/
├── input_files
│   └── 1000_Genomes_Project_chr16_unphased.vcf.gz
├── output_files
│   ├── 1000_Genomes_Project_chr16_unphased_SNV_genotypes.txt
│   └── prepare_genotype_file.log
├── prepare_genotype_file.py
└── README.md

```

3 directories, 5 files

---

## Overview

`prepare_genotype_file.py` is a Python utility to prepare SNV genotypes file from 
a vcf file for use in HaploGI runs.

---

## Requirements

- Python 3.x  
- No external dependencies (standard library only)

---

## Usage:

Run the script from the command line with the required inputs:

```bash
    python prepare_genotype_file.py --vcf <input.vcf or input.vcf.gz> 
                                   [--output_dir <directory>] 
```

---

## Arguments

| Argument       | Required | Description                                                                 |
|----------------|----------|-----------------------------------------------------------------------------|
| `--vcf`        | Yes      | Path to the input VCF file (`.vcf` or `.vcf.gz`) to create SNV genotypes file for HaploGI runs |
| `--output_dir` | No       | Directory where the output SNV genotypes file and log file will be saved (default: current directory) |

---

## Input Files

### 🔷 `1000_Genomes_Project_chr16_unphased.vcf.gz`

- Standard VCF file containing unphased genotypes
- Accepted input formats:
  - Uncompressed: `.vcf`
  - Compressed: `.vcf.gz`

---

## Output Files

### 🟢 `1000_Genomes_Project_chr16_unphased_SNV_genotypes.txt`

```bash
variant_position 302 302 303 303 306 306 307 307 402 402 403 403 404 404 406 406 407 407 408 408 410 410 411 411 412 412 414 414 415 415 416 416 501 501 502 502 503 503 504 504 505 505 506 506 507 507 508 508 509 509 510 510 511 511 512 512 513 513 514 514 515 515 516 516
16:10414 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
16:10638 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```

The SNV genotypes file contains variant genotype data for subjects with whole genome sequencing (WGS) data. Genotypes are encoded as:

- `0` = Missing  
- `1` = Reference allele (REF)  
- `2` = Alternative allele (ALT)  

#### File Structure

- **First column**:  
  Contains the SNV’s genomic position in the format `chromosome:position` (e.g., `16:10414`).

- **Remaining columns**:  
  Each subject is represented by **two consecutive columns**, one for each of their diploid genotype alleles.

- **Header row**:  
  Lists subject IDs. Each subject ID appears **twice**, corresponding to their two genotype alleles.

### 🟢 `prepare_genotype_file.log`
A log summary of the processing steps, including any warnings and errors encountered during execution.

---

## Example

Unzip the `input_files` directory before running the script. Use the included example inputs to run:

```bash
python prepare_genotype_file.py --vcf ./input_files/1000_Genomes_Project_chr16_unphased.vcf.gz \
							    --output_dir ./output_files_test
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
📧 nrafscience@gmail.com  
🔗 [GitHub Issues](https://github.com/RafPrograms/HaploGI/issues)

