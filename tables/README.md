# Table generation

This directory contains scripts used to generate or summarize tables reported in the manuscript.

## [`summarize_cooccurrence_network.py`](summarize_cooccurrence_network.py)

`summarize_cooccurrence_network.py` summarizes co-occurrence relationships among a predefined set of candidate variants and groups them into connected components.

The script was developed for analysis of candidate polymorphic markers surrounding the **EPCAM c.859-691_*1991del** deletion using variant co-occurrence information derived from the 54KJPN dataset of the Japanese Multi-Omics Reference Panel (jMorp).

## Requirements

* Python 3
* No external Python packages are required.

## External data dependency

Variant co-occurrence information was obtained from the Japanese Multi Omics Reference Panel (jMorp), provided by the Tohoku Medical Megabank Organization (ToMMo).

Dataset used in this study:

* Dataset: `Variation co-occurrence (54KJPN)`
* Genome build: `GRCh38`
* Release/date: `20230828`
* Date accessed: `2026-08-01`

The original jMorp dataset is not redistributed with this repository.

## Input

The input is a tab-delimited text file with a variable number of columns.

Each row represents one candidate variant:

```text
Position_of_variant    Position_of_cooccurred_variant1    Position_of_cooccurred_variant2    ...
```

Example:

```text
46644587    46651234    46658901    46672345
46821830    46830125
46955625    46960112    46967230
47089757
```

The first column contains the genomic position of the target variant. Subsequent columns contain genomic positions of variants showing co-occurrence with that target variant.

A header line is optional.

Empty fields and the values `.`, `NA`, and `N/A` are ignored.

## Usage

```bash
python summarize_cooccurrence_network.py cooccurrence_by_variant.tsv -o cooccurrence_network
```

## Output files

The script generates three tab-delimited output files.

### 1. `cooccurrence_network.summary.tsv`

Contains one row for each target variant.

Major columns include:

* `variant` — genomic position of the target variant
* `n_outgoing_all` — number of co-occurring variants listed in the input row
* `n_outgoing_target` — number of co-occurring variants that are also included in the predefined target set
* `n_incoming_target` — number of target variants that list the variant as a co-occurring partner
* `n_direct_target_partners` — total number of unique directly connected target variants
* `direction_status` — direction of the relationship in the input dataset
* `component_id` — connected-component identifier
* `component_size` — number of target variants in the connected component
* `direct_target_partners` — directly connected variants within the target set

`direction_status` is classified as:

* `BOTH` — relationships are observed in both outgoing and incoming directions
* `OUTGOING_ONLY` — co-occurring target variants are listed only in the row of the variant
* `INCOMING_ONLY` — the variant is listed only as a partner of another target variant
* `ISOLATED_WITHIN_TARGETS` — no direct co-occurrence relationship is observed with another variant in the predefined target set

**Note:** `ISOLATED_WITHIN_TARGETS` does not imply that the variant has no co-occurring variants in the complete 54KJPN dataset. It indicates only that no direct relationship was identified with another variant included in the input target set.

### 2. `cooccurrence_network.components.tsv`

Summarizes connected components among the target variants.

Columns include:

* `component_id`
* `component_size`
* `n_edges`
* `component_status`
* `members`

A component containing only one target variant is reported as `SINGLETON`.

### 3. `cooccurrence_network.edges.tsv`

Lists unique undirected edges between target variants.

Columns include:

* `variant1`
* `variant2`
* `direction`

The `direction` column retains information on whether the relationship was present in one or both directions in the original input.

## Network definition

Each target variant is treated as a node. A direct co-occurrence relationship between two target variants is treated as an edge.

Connected components are constructed only from variants present in the first column of the input file. Co-occurring variants outside the predefined target set are retained in the summary output but are not used to connect target variants indirectly.

Thus, if:

```text
A -> X
B -> X
```

and `X` is not included in the predefined target set, variants `A` and `B` are not considered directly connected by this script.

## Application in the EPCAM deletion analysis

The script was used to identify groups of polymorphic markers connected by direct co-occurrence relationships in the 54KJPN population dataset.

Variants belonging to the same connected component were considered potentially redundant when selecting representative markers. Variants classified as `ISOLATED_WITHIN_TARGETS` were subsequently examined directly in the original 54KJPN co-occurrence dataset to determine whether they had co-occurring variants outside the candidate set.

The final selection of candidate markers, including representative SNPs from LD clusters and isolated candidate variants, was performed separately from this script and is described in the Materials and Methods of the accompanying study.
