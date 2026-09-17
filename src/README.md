# `hapboundary4.py`

## Overview

`hapboundary4.py` estimates the shared disease-associated haplotype boundary between two phased samples around a specified genomic anchor.

The program was developed for pairwise haplotype analysis in which the disease-associated phased haplotype of a confirmed affected reference sample (`sample1`) is known.

The main concept is:

1. Fix the disease-associated haplotype of `sample1`.
2. Infer which phased haplotype of `sample2` matches `sample1` around the anchor.
3. Evaluate haplotype sharing in sliding windows extending in both directions from the anchor.
4. Report both strict and phase-switch-tolerant shared haplotype boundaries.

The switch-tolerant analysis allows the phased haplotype label of `sample2` to switch between B1 and B2 when persistent evidence supports a probable statistical phasing switch, while the disease-associated haplotype of `sample1` remains fixed.

---

## Requirements

`hapboundary4.py` uses only the Python standard library.

No external Python packages are required.

Required input files are:

- phased VCF or VCF.GZ
- genetic map
- family/sample metadata file

The VCF should contain phased diploid genotypes.

Only biallelic SNVs with genotype calls represented as:

```text
0|0
0|1
1|0
1|1
```

are used in the analysis.

---

## Basic usage

```bash
python hapboundary4.py \
  -i phased.vcf.gz \
  -m genetic_map.txt \
  -f family.txt \
  -a CHROM:POSITION \
  --sample1 SAMPLE_A \
  --sample2 SAMPLE_B \
  --anchor-allele 1 \
  -k 51 \
  --write-markers \
  -o output_prefix
```

Example:

```bash
python hapboundary4.py \
  -i phased.vcf.gz \
  -m genetic_map.txt \
  -f family.txt \
  -a 2:47611613 \
  --sample1 SAMPLE_A \
  --sample2 SAMPLE_B \
  --anchor-allele 1 \
  -k 51 \
  --write-markers \
  -o SAMPLE_A-SAMPLE_B
```

---

## Input files

### 1. Phased VCF

Specified with:

```text
-i / --input
```

The input may be either:

```text
VCF
VCF.GZ
```

The VCF must contain the samples specified by `--sample1` and `--sample2`.

The pathogenic variant or deletion itself does not need to be present in the VCF. The `--anchor` argument defines the genomic location around which the disease-associated haplotype is evaluated.

Only biallelic SNVs are used.

---

### 2. Genetic map

Specified with:

```text
-m / --map
```

The program uses the genetic map to interpolate genetic positions in centimorgans (cM).

Several common formats are recognized automatically, including formats containing combinations of:

```text
chromosome
physical position (bp)
genetic position (cM)
```

For example:

```text
2    46885797    67.123456
2    46886504    67.124001
2    46886604    67.124078
```

The program performs linear interpolation between adjacent map positions when an exact physical position is not present in the map.

Positions outside the range covered by the genetic map cannot be assigned an interpolated cM value.

---

### 3. Family/sample metadata

Specified with:

```text
-f / --family
```

The file must contain three fields:

```text
sample_id    family_id    affected
```

Header names are flexible. Common alternatives such as `IID`, `FID`, `status`, or `phenotype` are recognized.

Example:

```text
sample_id    family_id    affected
SAMPLE_A     FAMILY_1     1
SAMPLE_B     FAMILY_2     1
SAMPLE_C     FAMILY_2     0
```

Accepted values for the affected-status field include:

```text
1 / 0
affected / unaffected
case / control
yes / no
true / false
```

`sample1` must be marked as affected because it is treated as the confirmed disease-haplotype reference.

---

## Anchor and haplotype assignment

### Anchor position

The anchor is specified as:

```text
-a CHROM:POSITION
```

For example:

```text
-a 2:47611613
```

The anchor identifies the genomic location of the pathogenic variant or deletion.

The pathogenic variant itself may have been removed before phasing; the anchor is still used to identify the local haplotype surrounding its genomic position.

### Disease-associated haplotype of sample1

The user must specify which phased haplotype of `sample1` carries the disease-associated allele:

```text
--anchor-allele 1
```

or:

```text
--anchor-allele 2
```

The program then compares only the two candidate states consistent with this fixed `sample1` haplotype.

For example, if:

```text
--anchor-allele 1
```

the candidate states are:

```text
A1-B1
A1-B2
```

where A represents `sample1` and B represents `sample2`.

The corresponding haplotype of `sample2` is inferred from the marker-vote support in the anchor window.

---

## Marker voting

For each informative marker, the phased alleles of `sample1` and `sample2` are compared in four possible states:

```text
A1-B1
A1-B2
A2-B1
A2-B2
```

Marker votes are accumulated within a sliding window of `K` informative markers.

The default window size is:

```text
-k 51
```

`K` must be an odd integer of at least 3.

Two marker-selection modes are available.

### `original`

Default:

```text
--marker-mode original
```

This mode excludes identical homozygous sites but retains other informative marker configurations, including heterozygous-heterozygous and IBS0 sites.

### `hom-het`

```text
--marker-mode hom-het
```

This mode uses only homozygous-heterozygous marker combinations.

---

## Window classification

For each sliding window, support for the selected haplotype state is classified as:

```text
CONFIDENT_SHARED
POSSIBLE_SHARED
NON_SHARED
```

Default thresholds are:

```text
CONFIDENT_SHARED: selected-state fraction >= 2/3
POSSIBLE_SHARED:  selected-state fraction >= 1/2
NON_SHARED:       below the possible threshold
```

These thresholds can be modified using:

```text
--confident-fraction
--possible-fraction
```

---

## Boundary definitions

The program reports four haplotype-boundary estimates.

### `strict_confident`

A contiguous region around the anchor classified only as:

```text
CONFIDENT_SHARED
```

using the initially selected sample1-sample2 haplotype state without allowing phase-label changes.

### `strict_possible`

A contiguous region around the anchor classified as either:

```text
CONFIDENT_SHARED
```

or:

```text
POSSIBLE_SHARED
```

without allowing phase-label changes.

### `switch_tolerant_confident`

Equivalent to the confident boundary, but permits B1/B2 label switching in `sample2` when the alternative state is persistently supported.

### `switch_tolerant_possible`

Equivalent to the possible boundary, while allowing accepted B1/B2 label switches in `sample2`.

The `sample1` disease-associated haplotype remains fixed throughout the switch-tolerant analysis.

---

## Phase-switch tolerance

Statistical phasing may cause the B1/B2 labels of `sample2` to switch across a genomic region.

`hapboundary4.py` can recognize such transitions when support for the alternative `sample2` haplotype is sufficiently strong and persistent.

Default parameters are:

```text
--switch-min-fraction  0.6667
--switch-min-margin    5
--switch-min-windows   5
```

A switch is accepted only when the alternative state:

1. reaches the required support fraction,
2. exceeds the current state by the required vote margin, and
3. satisfies these criteria across the specified number of consecutive windows.

Left and right directions are evaluated independently from the anchor.

---

## Anchor-assignment parameters

The assignment of the matching `sample2` haplotype at the anchor is controlled by:

```text
--anchor-min-fraction
--anchor-min-margin
```

Defaults:

```text
--anchor-min-fraction 0.5
--anchor-min-margin   1
```

If there is insufficient support or the two candidate haplotypes cannot be distinguished, the program terminates with an anchor-assignment error rather than assigning an ambiguous haplotype.

---

## Output files

For an output prefix:

```text
-o SAMPLE_A-SAMPLE_B
```

the program generates the following files.

### `<prefix>.summary.tsv`

Main summary output containing:

- sample and family information
- affected status
- anchor position
- inferred `sample2` anchor haplotype
- selected haplotype state
- anchor-window statistics
- analysis parameters
- accepted phase switches
- strict confident boundary
- strict possible boundary
- switch-tolerant confident boundary
- switch-tolerant possible boundary

Both physical length (bp) and genetic length (cM) are reported where available.

---

### `<prefix>.windows.tsv`

Window-level output containing:

- genomic and genetic coordinates of each sliding window
- vote counts for all four haplotype states
- strict selected state and classification
- switch-tolerant selected state and classification
- detected switch events
- IBS0 counts

This file can be used to inspect how haplotype support changes across the analyzed region.

---

### `<prefix>.switches.tsv`

Lists accepted phase-label switch events in `sample2`.

Reported fields include:

- direction from anchor
- genomic position
- genetic position
- original state
- alternative state
- vote counts
- vote margin
- number of supporting windows

---

### `<prefix>.markers.tsv`

Generated only when:

```text
--write-markers
```

is specified.

This file reports marker-level information including:

- position
- interpolated cM
- phased genotypes of both samples
- marker type
- votes for A1-B1, A1-B2, A2-B1, and A2-B2
- selected-state vote
- IBS0 status

---

## Main command-line options

| Option | Description | Default |
|---|---|---|
| `-i`, `--input` | Phased VCF or VCF.GZ | required |
| `-m`, `--map` | Genetic map | required |
| `-f`, `--family` | Family/sample metadata | required |
| `-a`, `--anchor` | Anchor coordinate (`CHROM:POS`) | required |
| `--sample1` | Confirmed affected reference sample | required |
| `--sample2` | Sample compared with sample1 | required |
| `--anchor-allele` | Disease-carrying phased haplotype of sample1 (1 or 2) | required |
| `-k` | Number of informative markers per sliding window | 51 |
| `--confident-fraction` | Threshold for `CONFIDENT_SHARED` | 2/3 |
| `--possible-fraction` | Threshold for `POSSIBLE_SHARED` | 0.5 |
| `--anchor-min-fraction` | Minimum support for anchor assignment | 0.5 |
| `--anchor-min-margin` | Minimum vote difference at anchor | 1 |
| `--switch-min-fraction` | Minimum support for sample2 phase switch | 2/3 |
| `--switch-min-margin` | Minimum vote margin for phase switch | 5 |
| `--switch-min-windows` | Consecutive qualifying windows required | 5 |
| `--marker-mode` | Marker-selection mode (`original` or `hom-het`) | `original` |
| `--write-markers` | Write marker-level output | off |
| `-o`, `--output-prefix` | Output prefix | required |

The complete command-line help can also be displayed using:

```bash
python hapboundary4.py --help
```

---

## Interpretation

The strict boundaries assume that the phased B1/B2 labels of `sample2` remain consistent across the entire analyzed interval.

The switch-tolerant boundaries are intended to distinguish a true loss of shared haplotype from an apparent boundary caused by a probable phase-label switch in `sample2`.

The choice of boundary metric used for downstream analyses should be specified in the associated manuscript and analysis documentation.

---

## Relationship to the full workflow

The upstream preparation of phased VCF files and the commands used in the study are described in:

[`../docs/analysis_commands.md`](../docs/analysis_commands.md)

The main analysis script is:

[`hapboundary4.py`](hapboundary4.py)

Figure-generation scripts are documented in:

[`../figures/README.md`](../figures/README.md)