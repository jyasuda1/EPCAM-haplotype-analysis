# Figure generation

This directory contains the scripts used to generate the figures reported in the manuscript.

## External data dependencies


Generation of some figures requires reference data obtained from the Japanese Multi Omics Reference Panel (jMorp), provided by the Tohoku Medical Megabank Organization (ToMMo).

The original jMorp data are not included in this repository. Users should obtain the corresponding datasets directly from jMorp and follow the applicable terms and conditions of use.

### `plot_recombination_rate.R`

This script uses recombination-rate data obtained from jMorp.

Dataset used in this study:

* Source: Japanese Multi Omics Reference Panel (jMorp)
* Provider: Tohoku Medical Megabank Organization (ToMMo)
* Dataset: `Genetic Map`
* Genome build: `GRCh37`
* Release/date: `20210907`
* Date accessed: `2026-08-01`

The downloaded data were processed into the input format described below before plotting.

### [`make_candidate_r2_matrix.py`](../scripts/make_candidate_r2_matrix.py)

This script uses variant co-occurrence data obtained from jMorp to generate the pairwise r² matrix used for Figure 4B.

Dataset used in this study:

* Source: Japanese Multi Omics Reference Panel (jMorp)
* Provider: Tohoku Medical Megabank Organization (ToMMo)
* Dataset: `Variation co-occurrence (54KJPN)`
* Genome build: `GRCh38`
* Release/date: `20230828`
* Date accessed: `2026-08-01`

The original jMorp dataset is not redistributed with this repository.

---

## Figure 3A

**Script**

`plot_haplotype_core_flank.R`

**Purpose**

Generates the haplotype distribution plot showing the shared core haplotype region and flanking haplotype regions for the analyzed individuals.

**Input**

A tab-delimited file containing the genomic or genetic boundaries of the shared haplotype regions and plotting information for each sample.

Expected columns:

```text
sample_name    start_position    end_position    core_color    flanking_color
```

Additional details regarding the input format and plotting parameters are provided in the script.

**Output**

Figure 3A as PNG and PDF files.

**Requirements**

R with the packages specified in the script.

---

## Figure 3B

**Script**

`plot_recombination_rate.R`

**Purpose**

Generates the recombination-rate plot surrounding the haplotype core region.

**Input**

A genetic map file containing physical genomic positions and local recombination rates.

Example format:

```text
position    recombination_rate
46885797    0.773550266
46886504    0.746729130
46886604    0.727462649
```

**Output**

Figure 3B as PNG and PDF files.

**Requirements**

R with the packages specified in the script.

---

## Figure 4B

**Scripts**

[`make_candidate_r2_matrix.py`](../scripts/make_candidate_r2_matrix.py)
[`plot_r2_heatmap_ldblocks_intercluster.R`](plot_r2_heatmap_ldblocks_intercluster.R)

**Purpose**

Generates the pairwise linkage-disequilibrium matrix and the corresponding heatmap used in Figure 4B.

**Workflow**

```text
jMorp Variation co-occurrence data
    |
    v
make_candidate_r2_matrix.py
    |
    v
pairwise r² matrix
    |
    v
plot_r2_heatmap_ldblocks_intercluster.R
    |
    v
Figure 4B
```

**Input**

`make_candidate_r2_matrix.py` takes the jMorp 54KJPN variant co-occurrence dataset and a list of candidate SNP positions as input and generates a pairwise r² matrix. The resulting matrix is used as input for `plot_r2_heatmap_ldblocks_intercluster.R`.

Details of the required input format and parameters are provided in the respective scripts.


**Output**

Figure 4B as PNG and PDF files.

**Requirements**

Python and R with the packages specified in the respective scripts.

---

## Reproducibility

The scripts in this directory correspond to the final analyses used to generate the figures reported in the manuscript.

For details of the upstream analyses, see the [main README](../README.md), [`hapboundary4.py`](../src/hapboundary4.py), and [`analysis_commands.md`](../docs/analysis_commands.md).
