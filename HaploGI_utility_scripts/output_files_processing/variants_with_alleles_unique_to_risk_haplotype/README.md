# `risk_alleles_variants.py`

![Python](https://img.shields.io/badge/python-3.x-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3.0-blue)

🔗 [`risk_alleles_variants.py`](./risk_alleles_variants.py)

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
variants_with_alleles_unique_to_risk_haplotype/
├── input_files
│   ├── 1000_Genomes_Project_chr16_unphased.vcf.gz
│   ├── chr16_SNV_genotypes_B38.txt
│   ├── core_set_of_cases.txt
│   └── haplotype_sharing_controls.txt
├── output_files
│   ├── risk_alleles_variants.log
│   └── risk_alleles_variants.txt
├── README.md
└── risk_alleles_variants.py

```

3 directories, 8 files

Folders containing input and output files are zipped to save space.

---

## Overview

`risk_alleles_variants.py` is a Python utility for identifying variants within a risk haplotype that carry alleles unique to that haplotype—referred to as *risk alleles*. 

The script analyzes genotype data to find variants that meet the following criteria:
- All individuals in the core set of cases are **heterozygous** at the variant site.
- All other individuals (cases and controls **not** carrying the risk haplotype) are **homozygous** for either the reference or alternate allele.
- Controls who share a significant portion of the risk haplotype are excluded from this analysis to avoid confounding.

This tool helps isolate potentially causal variants uniquely associated with the risk haplotype for downstream interpretation or validation.

---

## Requirements
 
- Python 3.6+
- pandas

---

## Usage

Run the script from the command line with the required inputs:

```bash
 python risk_alleles_variants.py --start_pos <int> \
                                 --end_pos <int> \
                                 --genotypes <genotype_file> \
                                 --set <core_cases_file> \
                                 --vcf <vcf_file(.vcf or .vcf.gz)> \
                                [--controls <shared_haplotype_controls_file>] \
                                [--output_dir <directory>]
```
 
---

## Arguments

| Argument        | Required | Description                                                                                             |
|-----------------|----------|---------------------------------------------------------------------------------------------------------|
| `--start_pos`   | Yes      | Starting base pair position at the left boundary of the risk haplotype region.                          |
| `--end_pos`     | Yes      | Ending base pair position at the right boundary of the risk haplotype region.                           |
| `--genotypes`   | Yes      | Path to the HaploGI genotype file (e.g., `SNV_genotypes.txt`).                                          |
| `--set`         | Yes      | Path to the core set of cases file (a single line listing subject IDs who share the risk haplotype).    |
| `--vcf`         | Yes      | Path to the original VCF file used to generate the genotype file (`.vcf` or `.vcf.gz`).                 |
| `--controls`    | No       | Path to a file listing control subjects who also share the risk haplotype (optional).                   |
| `--output_dir`  | No       | Directory where output files will be saved (default: current directory).                                |

---

## Input Files

### 🔷 `1000_Genomes_Project_chr16_unphased.vcf.gz`
- Standard VCF file containing unphased genotypes.
- Used as input to generate the SNV genotypes file for HaploGI runs.

### 🔷 `chr16_SNV_genotypes_B38.txt`
```bash
variant_position 302 302 303 303 306 306 307 307 402 402 403 403 404 404 406 406 407 407 408 408 410 410 411 411 412 412 414 414 415 415 416 416 501 501 502 502 503 503 504 504 505 505 506 506 507 507 508 508 509 509 510 510 511 511 512 512 513 513 514 514 515 515 516 516
16:10414 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
16:10638 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```

The SNV genotype file contains variant genotype data for subjects with whole genome sequencing (WGS) data. Genotypes are encoded as:

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


### 🔷 `core_set_of_cases.txt`
```bash
302 306 403 408 411 501 504 506 511 512 513 516
```

This file contains a list of **case subject IDs**, separated by spaces, all on a **single line**.

- **No header row**

HaploGI uses this set of cases to **check for the existence of haplotype sharing** among them.

### 🔷 `haplotype_sharing_controls.txt`
```bash
404 503
```
- This file contains a list of **control subject IDs** who share a significant portion of a risk haplotype.
- IDs are space-separated and appear on a single line.

- **No header row**

---

## Output Files

### 🟢 `risk_alleles_variants.txt`
```bash
Position	Risk_Allele(0=REF;1=ALT)	ID	            REF	ALT
1895932	    1	                    16:1895932:A:G	A	G
1957759	    1	                    16:1957759:G:C	G	C
1973742	    1	                    16:1973742:C:T	C	T
2036534	    1	                    16:2036534:A:G	A	G
```

- Position: Base pair position of the variant with an allele unique to the risk haplotype

- Risk_Allele (0=REF; 1=ALT): Indicates which allele (REF or ALT) is unique to the risk haplotype

- ID, REF, ALT: Metadata from the original VCF file used to generate the SNV genotype file for HaploGI runs

### 🟢 `risk_alleles_variants.log`
A log summary of the processing steps, including any warnings and errors encountered during execution.

---

## Example

Unzip the `input_files` directory before running the script. Use the included example inputs to run:

```bash
python risk_alleles_variants.py --start_pos 1831445 \
                 --end_pos 8321356 \
                 --genotypes ./input_files/chr16_SNV_genotypes_B38.txt \
                 --set ./input_files/core_set_of_cases.txt \
                 --vcf ./input_files/1000_Genomes_Project_chr16_unphased.vcf.gz \
                 --controls ./input_files/haplotype_sharing_controls.txt \
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

See the [LICENSE](./LICENSE) file in this repository for complete details.

---

## Support

For questions, bug reports, or suggestions, please contact:  
📧 your.email@example.com  
🔗 [GitHub Issues](https://github.com/yourusername/yourrepo/issues)

