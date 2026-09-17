# Haplotype analysis of the 4.13-kb EPCAM deletion in Japanese Lynch syndrome families

## Overview

This repository contains custom scripts and computational procedures used for haplotype analysis of the 4.13-kb EPCAM deletion, a founder mutation identified in Japanese Lynch syndrome families.

The workflow includes SNP-array preprocessing, haplotype phasing, estimation of shared haplotype boundaries, and downstream analyses based on variant co-occurrence data from the Japanese population.

The main custom program, `hapboundary4.py`, was developed to estimate shared haplotype boundaries surrounding the EPCAM deletion locus.

## Associated publication

This repository accompanies the following manuscript currently in preparation:

Haplotype analysis of the 4.13-kb EPCAM deletion in Japanese Lynch syndrome families

Yasuda et al.
Manuscript in preparation.

## Repository structure

```text
EPCAM_project/
├── README.md
├── CITATION.cff
├── src/
│   ├── README.md
│   └── hapboundary4.py
├── scripts/
│   └── shapeit_haps_to_vcf.py
├── figures/
│   ├── README.md
│   └── ...
├── tables/
│   ├── README.md
│   └── ...
└── docs/
    └── analysis_commands.md
```

* `src/` - main custom analysis code and documentation for `hapboundary4.py`.
* `scripts/` - utility scripts used for file conversion and preprocessing.
* `figures/` - scripts and documentation used to generate manuscript figures.
* `tables/` - scripts and documentation used to generate or summarize manuscript tables.
* `docs/` - detailed analysis commands and the computational workflow used in the study.


## Analysis workflow

```text
Japonica array (a SNP array for Japanese population genetics)
        |
        v
Quality control / annotation correction
        |
        v
SHAPEIT phasing
        |
        v
Conversion to phased VCF
        |
        v
REF/ALT orientation correction
        |
        v
hapboundary4.py
        |
        +--> shared haplotype boundaries
        |
        +--> haplotype figures
        |
        +--> core haplotype SNP selection
                  |
                  v
           GRCh37-to-GRCh38 LiftOver
                  |
                  v
           jMorp LD / co-occurrence analysis
```

Detailed commands are provided in
[`docs/analysis_commands.md`](docs/analysis_commands.md).

## Main analysis program

### `hapboundary4.py`

`hapboundary4.py` estimates shared haplotype boundaries between phased haplotypes around a specified anchor variant or locus.

Main inputs:

* phased VCF
* genetic map
* family/sample information
* anchor position
* disease-associated phased haplotype (1 or 2) of the reference sample

The program reports several definitions of shared haplotype boundaries, including strict and switch-tolerant estimates.

For the analyses reported in the manuscript, the `strict_possible` boundary was used as the representative shared haplotype interval, whereas the `switch_tolerant_possible` boundary was used to identify candidate proxy markers for the EPCAM deletion.

Detailed usage, algorithmic definitions, input requirements, and output formats are provided in:

[`src/README.md`](src/README.md)

The methodological rationale for the hapboundary algorithm and the study-specific parameter settings are described in the Supplementary Methods of the associated manuscript.


## External data sources

Reference genetic map and variant co-occurrence data were obtained from the Japanese Multi Omics Reference Panel (jMorp), provided by the Tohoku Medical Megabank Organization (ToMMo).

The original jMorp datasets are not redistributed in this repository.

Dataset-specific information is provided in the corresponding documentation in the `figures/`, `tables/`, and `docs/` directories.

Additional reference resources used for variant annotation and coordinate conversion are documented in `docs/analysis_commands.md`.

## Software requirements

The computational workflow used the following software:

- PLINK v1.9
- SHAPEIT v2.r904
- bcftools v1.23.1
- Python 3
- R v4.5.2
- ANNOVAR (2020-06-07 release)
- HTSlib (bgzip/tabix) v1.23.1

## Reproducibility

The commands used for SNP-array preprocessing, phasing, VCF generation, haplotype analysis, and downstream LD analyses are documented in: [`docs/analysis_commands.md`](docs/analysis_commands.md)

## Data availability

Individual-level genotype and clinical data are not included in this repository.

Reference datasets used in the analyses should be obtained directly from their original providers and are subject to the corresponding terms of use.

## Code availability

Custom Python and R scripts used in the manuscript are provided in this repository.

## Citation

If you use the code in this repository, please cite the associated publication.

Citation information is also provided in `CITATION.cff`.
