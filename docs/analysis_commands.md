# Analysis commands

This document summarizes the commands and processing steps used for SNP-array preprocessing, phasing, haplotype-boundary analysis, and downstream linkage-disequilibrium analyses.

The commands are provided as a record of the computational workflow used in the study. File names correspond to those used during the original analysis.

## Environment-specific paths

The following paths should be modified according to the local computing environment.

```bash
REF_FASTA=/path/to/hs37d5.fa
BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
ANNOVAR_DB=/path/to/annovar/humandb
```

---

# Step 1. Preparation and pruning of Japonica array data

## 1.1 Remove poor-quality SNP probes

A list of poor probes was generated from the SNP summary table.

```bash
awk '{FS=OFS="\t"}$21 ~/^O/||$21 ~/^C/ {print $1}' \
  SNP_summary_table_v2.txt > poorprove.txt
```

The corresponding probes were excluded using PLINK.

```bash
plink \
  --file genotyping_result_v2 \
  --exclude poorprove.txt \
  --make-bed \
  --out Lynch1
```

If PLINK returns an error caused by extra chromosome identifiers, the `--aec` option may be added where appropriate.

## 1.2 Update sample and pedigree information

The following files were used to update individual information:

* `updateID.txt`
* `updateParent.txt`
* `updateSex.txt`
* `updatephenotype.txt`

The file formats follow the PLINK individual-data update specifications.

### Update sample IDs

```bash
plink \
  --bfile Lynch1 \
  --update-ids updateID.txt \
  --allow-no-sex \
  --make-bed \
  --recode \
  --out Lynch2
```

### Update parental information

```bash
plink \
  --bfile Lynch2 \
  --update-parents updateParent.txt \
  --allow-no-sex \
  --make-bed \
  --recode \
  --out Lynch3
```

### Update sex

```bash
plink \
  --bfile Lynch3 \
  --update-sex updateSex.txt \
  --recode \
  --make-bed \
  --out Lynch4
```

### Update phenotype

```bash
plink \
  --bfile Lynch4 \
  --pheno updatephenotype.txt \
  --recode \
  --make-bed \
  --out Lynch5
```

## 1.3 Quality-control analyses

### Check sex

```bash
plink \
  --bfile Lynch5 \
  --check-sex
```

### Remove variants with MAF < 0.02

```bash
plink \
  --bfile Lynch5 \
  --maf 0.02 \
  --make-bed \
  --out Lynch6
```

### Check cryptic relatedness

```bash
plink \
  --bfile Lynch6 \
  --genome \
  --min 0.015625 \
  --out p-hatLynch6
```

## 1.4 Principal-component analysis

### Missingness and heterozygosity

```bash
plink \
  --bfile Lynch6 \
  --missing \
  --out qc_missingLynch6

plink \
  --bfile Lynch6 \
  --het \
  --out qc_hetLynch6
```

### Select variants for PCA

```bash
plink \
  --bfile Lynch6 \
  --autosome \
  --maf 0.05 \
  --geno 0.02 \
  --mind 0.02 \
  --make-bed \
  --out Lynch6_qc
```

LD pruning was performed as follows:

```bash
plink \
  --bfile Lynch6_qc \
  --indep-pairwise 50 5 0.2 \
  --out Lynch6_pruned
```

PCA was then performed using the pruned variants.

```bash
plink \
  --bfile Lynch6_qc \
  --extract Lynch6_pruned.prune.in \
  --pca 10 \
  --out PCA
```

### PCA plot

```r
library(ggplot2)
library(ggrepel)

pca <- read.table(
    "PCA.eigenvec",
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(pca)[1:4] <- c("FID", "IID", "PC1", "PC2")

pdf("wholememberPCA.pdf")

ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(size = 3) +
    geom_text_repel(
        aes(label = IID),
        size = 3.5,
        max.overlaps = Inf
    ) +
    theme_classic() +
    labs(
        x = "PC1",
        y = "PC2"
    )

dev.off()
```

## 1.5 PCA using one representative individual per family

The representative individuals were listed in `PCA_representative.txt`.

Example:

```text
F1      W1423
F2      W1742
W308    W308
W732    W732
W1349   W1349
W1324   W1324
W695    W695
W945    W945
GT650   GT650
W1663   W1663
GT441   GT441
GT399   GT399
GT484   GT484
```

PCA was performed as follows:

```bash
plink \
  --bfile Lynch6_qc \
  --keep PCA_representative.txt \
  --extract Lynch6_pruned.prune.in \
  --pca 10 \
  --out PCA_select
```

The PCA plot was generated in R.

```r
library(ggplot2)
library(ggrepel)

pca <- read.table(
    "PCA_select.eigenvec",
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(pca)[1:4] <- c("FID", "IID", "PC1", "PC2")

pdf("selectedPCA.pdf")

ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(size = 3) +
    geom_text_repel(
        aes(label = IID),
        size = 3.5,
        max.overlaps = Inf
    ) +
    theme_classic() +
    labs(
        x = "PC1",
        y = "PC2"
    )

dev.off()
```

A second representative set was analyzed using `PCA_representative2.txt`.

```bash
plink \
  --bfile Lynch6_qc \
  --keep PCA_representative2.txt \
  --extract Lynch6_pruned.prune.in \
  --pca 10 \
  --out PCA_select2
```

```r
library(ggplot2)
library(ggrepel)

pca <- read.table(
    "PCA_select2.eigenvec",
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(pca)[1:4] <- c("FID", "IID", "PC1", "PC2")

pdf("selectedPCA2.pdf")

ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(size = 3) +
    geom_text_repel(
        aes(label = IID),
        size = 3.5,
        max.overlaps = Inf
    ) +
    theme_classic() +
    labs(
        x = "PC1",
        y = "PC2"
    )

dev.off()
```

## 1.6 Preparation of chromosome 2 data

When the analysis was transferred to another computing environment, the dataset prefix was renamed from `Lynch6` to `Lynch7`.

No genotype-processing step was performed as part of this renaming.

### Remove duplicate variants

```bash
plink \
  --bfile Lynch7 \
  --list-duplicate-vars \
  --out dup7

plink \
  --bfile Lynch7 \
  --exclude dup7.dupvar \
  --recode \
  --make-bed \
  --out Lynch8
```

### Select chromosome 2

```bash
plink \
  --bfile Lynch8 \
  --chr 2 \
  --geno 0.2 \
  --recode \
  --make-bed \
  --out Lynch8-chr2
```

### Remove variants with MAF < 0.02

```bash
plink \
  --bfile Lynch8-chr2 \
  --maf 0.02 \
  --make-bed \
  --out Lynch9-chr2
```

---

# Step 2. Phasing

The EPCAM deletion interval was excluded before phasing.

```bash
plink \
  --bfile Lynch9-chr2 \
  --e 2:47611163-47615743 \
  --make-bed \
  --out Lynch9-chr2delremove
```

Parental information was updated.

```bash
plink \
  --bfile Lynch9-chr2delremove \
  --update-parents updateParent.txt \
  --make-bed \
  --out Lynch9-chr2delremoveUF
```

Variants and samples were further filtered before phasing.

```bash
plink \
  --bfile Lynch9-chr2delremoveUF \
  --me 1 0 \
  --geno 0 \
  --mendel-duos \
  --mind 0.1 \
  --snps-only \
  --make-bed \
  --out Lynch9-chr2delremove_clean
```

## Genetic map

The genetic map was downloaded from the Japanese Multi Omics Reference Panel (jMorp):

https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-genetic_map-20210907/files/20210907tommo_genetic_map_300haploid.zip

The chromosome 2 genetic map file,

```text
genetic_map_tommo_GRCh37_chr2_bp5.txt
```

was manually converted to the map format required by SHAPEIT.

The resulting file was:

```text
ToMMochr2shapeit.map
```

## SHAPEIT phasing

```bash
shapeit \
  --input-bed Lynch9-chr2delremove_clean \
  --input-map ToMMochr2shapeit.map \
  --window 0.1 \
  --states 50 \
  --duohmm \
  --output-max \
  Lynch9-chr2delremove.phased.haps \
  autosome.chr2.imputed.plink.sample
```

---

# Step 3. Generation of a phased VCF file

The custom utility script [`shapeit_haps_to_vcf.py`](../scripts/shapeit_haps_to_vcf.py) was used to convert the SHAPEIT `.haps` output to VCF format.

```bash
python3 shapeit_haps_to_vcf.py \
  Lynch9-chr2delremove.phased.haps \
  Lynch9-chr2rmdel.phased.haps.noheader.vcf
```

A VCF header was manually generated and stored as:

```text
Lynch9vcfheader.txt
```

The header and variant records were concatenated.

```bash
cat \
  Lynch9vcfheader.txt \
  Lynch9-chr2rmdel.phased.haps.noheader.vcf \
  > Lynch9-chr2rmdel.vcf
```

The VCF was compressed and indexed using bcftools.

```bash
bcftools view \
  Lynch9-chr2rmdel.vcf \
  -Oz \
  -o Lynch9-chr2rmdel.vcf.gz

bcftools index \
  -t Lynch9-chr2rmdel.vcf.gz
```

---

# Step 4. Correction of REF/ALT strand-orientation inconsistencies

The bcftools plugin directory was exported before running `+fixref`.

```bash
export BCFTOOLS_PLUGINS
```

Reference consistency was checked against the GRCh37 reference FASTA.

```bash
bcftools +fixref \
  Lynch9-chr2rmdel.vcf.gz \
  -- \
  -f "$REF_FASTA"
```

The dbSNP GRCh37 VCF files used for reference annotation were obtained from NCBI.

```bash
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi
```

REF/ALT orientation was corrected using `bcftools +fixref`.

```bash
bcftools +fixref \
  Lynch9-chr2rmdel.vcf.gz \
  -Oz \
  -o Lynch9-chr2rmdelfixproto.vcf.gz \
  -- \
  -f "$REF_FASTA" \
  -m flip \
  -d
```

The corrected VCF was indexed, sorted, and re-indexed.

```bash
bcftools index \
  -t Lynch9-chr2rmdelfixproto.vcf.gz

bcftools sort \
  Lynch9-chr2rmdelfixproto.vcf.gz \
  -Oz \
  -o Lynch9-chr2rmdelfix.vcf.gz

bcftools index \
  -t Lynch9-chr2rmdelfix.vcf.gz
```

ANNOVAR was used with an uncompressed VCF input.

```bash
gunzip -c \
  Lynch9-chr2rmdelfix.vcf.gz \
  > Lynch9-chr2rmdelfix.vcf
```

Variant annotation was performed using RefGene and 54KJPN.

```bash
perl table_annovar.pl \
  Lynch9-chr2rmdelfix.vcf \
  "$ANNOVAR_DB" \
  --buildver hg19 \
  --out Lynch9-chr2rmdelfix \
  --protocol refGene,54KJPN \
  --operation g,f \
  --remove \
  --otherinfo \
  --vcfinput \
  --nastring .
```

---

# Step 5. Haplotype-boundary analysis

Haplotype boundaries were estimated using the custom Python program [`hapboundary4.py`](../src/hapboundary4.py).

The analysis was performed for all pairwise combinations of the study samples and both anchor alleles.

```bash
for i in \
  F1_W16 F1_W1423 F1_W1508 F1_W1509 F1_W1510 F1_W1511 \
  F1_W1521 F1_W1537 F1_W1538 F1_W1539 F1_W1551 F1_W1595 \
  F1_W1596 F2_W1653 F2_W1742 W308_W308 W732_W732 W1324_W1324 \
  W695_W695 W945_W945 W1349_W1349 GT650_GT650 W1663_W1663 \
  GT441_GT441 GT399_GT399 MLH_GT484 MLH_W1590

do

  for j in \
    F1_W16 F1_W1423 F1_W1508 F1_W1509 F1_W1510 F1_W1511 \
    F1_W1521 F1_W1537 F1_W1538 F1_W1539 F1_W1551 F1_W1595 \
    F1_W1596 F2_W1653 F2_W1742 W308_W308 W732_W732 W1324_W1324 \
    W695_W695 W945_W945 W1349_W1349 GT650_GT650 W1663_W1663 \
    GT441_GT441 GT399_GT399 MLH_GT484 MLH_W1590

  do

    for k in 1 2

    do

      python hapboundary4.py \
        -i Lynch9-chr2rmdelfix.vcf.gz \
        -m chr2geneticmapTommo3.txt \
        -f family.txt \
        -a 2:47611613 \
        --sample1 "${i}" \
        --sample2 "${j}" \
        -k 51 \
        --anchor-allele "${k}" \
        --write-markers \
        -o "${i}-${j}_${k}"

      wait

    done

  done

done
```

---

# Step 6. Variant co-occurrence analysis

The jMorp 54KJPN variant co-occurrence data for chromosome 2 were downloaded as:

```text
tommo-54kjpn-20230828-GRCh38-autosome-chr2-plink-r2.tsv.gz
```

The file `hg38pos157.txt` contains 157 SNPs selected from the maximal `switch_tolerant_possible` core haplotype interval identified using GRCh37 coordinates.

Their genomic positions were converted from GRCh37 to GRCh38 using UCSC LiftOver, producing:

```text
hg38pos157.txt
```

The custom Python script [`make_candidate_r2_matrix.py`](../scripts/make_candidate_r2_matrix.py) was used to generate the pairwise r族 matrix.

```bash
python make_candidate_r2_matrix.py \
  tommo-54kjpn-20230828-GRCh38-autosome-chr2-plink-r2.tsv.gz \
  hg38pos157.txt \
  -o pos157
```

## Identify variants without connections in the initial r族 matrix

```bash
awk '$4 == "NA" {print $1}' \
  pos157.r2_summary.tsv \
  > pos157.r2NA.tsv
```

Incoming connections to these variants were then checked directly against the jMorp co-occurrence dataset.

```bash
awk -F'\t' '
BEGIN {
    OFS = "\t"
}

NR == FNR {
    target[$1] = 1
    next
}

$9 in target {
    count[$9]++
}

END {
    print "Position", "n_incoming"

    for (p in target)
        print p, count[p] + 0
}
' pos157.r2NA.tsv \
  <(zcat tommo-54kjpn-20230828-GRCh38-autosome-chr2-plink-r2.tsv.gz) \
  | sort -k1,1n \
  > isolated_incoming.count157nolink.tsv
```

## Generate the LD heatmap

`SNPlist.txt` is a one-column text file containing the chromosome 2 positions of the 157 SNPs that were positive in all affected alleles within the core haplotype, using GRCh38 coordinates.

The coordinate

```text
47384474
```

corresponds to the EPCAM deletion position in GRCh38.

The heatmap was generated using [`plot_r2_heatmap_ldblocks_intercluster.R`](../figures/plot_r2_heatmap_ldblocks_intercluster.R).

```bash
Rscript plot_r2_heatmap_ldblocks_intercluster.R \
  pos157.r2_matrix.tsv \
  SNPlist.txt \
  47384474
```

---

# Step 7. Annotation of 54KJPN variants in the LD block

`hg38_54KJPN.txt` is an ANNOVAR database file generated previously from the 54KJPN dataset.

The database file was sorted and compressed.

```bash
(
    grep '^#' hg38_54KJPN.txt
    grep -v '^#' hg38_54KJPN.txt | sort -k1,1 -k2,2n -k3,3n
) | bgzip -c > hg38_54KJPN.txt.gz
```

A tabix index was generated.

```bash
tabix \
  -s 1 \
  -b 2 \
  -e 3 \
  -c '#' \
  hg38_54KJPN.txt.gz
```

Variants within the EPCAM LD interval were extracted.

```bash
tabix \
  hg38_54KJPN.txt.gz \
  chr2:47346929-47453034 \
  | awk 'BEGIN {FS=OFS="\t"} $4 != "-" && $5 != "-" && length($4) == 1 && $7 == "PASS" {print $0}' \
  > hg38_54KJPNEPCAMLD.txt
```

The extracted variants were annotated using ANNOVAR.

```bash
perl table_annovar.pl \
  hg38_54KJPNEPCAMLD.txt \
  "$ANNOVAR_DB" \
  --buildver hg38 \
  --out Lynch9-chr2rmdelfix \
  --protocol refGene,clinvar_20250706 \
  --operation g,f \
  --remove \
  --otherinfo \
  --nastring .
```

## LiftOver to GRCh37

The LD interval was converted from GRCh38 to GRCh37 using the UCSC Genome Browser LiftOver tool.

```text
GRCh38:
chr2:47346929-47453034

GRCh37:
chr2:47574067-47680173
```
