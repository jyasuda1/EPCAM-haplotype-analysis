#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_r2_heatmap_ldblocks.R R2_matrix.tsv SNPlist.txt anchor_or_anchorfile")
}

matrix_file <- args[1]
snp_file <- args[2]
anchor_arg <- args[3]

R2_THRESHOLD <- 0.50
MIN_BLOCK_SNPS <- 3
LABEL_EVERY <- 1

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2: install.packages('ggplot2')")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Please install scales: install.packages('scales')")
}
library(ggplot2)

read_anchor <- function(x) {
  if (file.exists(x)) {
    z <- read.table(x, header=FALSE, stringsAsFactors=FALSE)
    a <- suppressWarnings(as.numeric(z[[1]]))
    if (any(is.na(a))) stop("Anchor file contains non-numeric values.")
    return(a)
  }
  a <- suppressWarnings(as.numeric(x))
  if (is.na(a)) stop("Third argument must be an integer coordinate or a readable file.")
  a
}

coord_to_index <- function(coord, positions) {
  if (coord < min(positions) || coord > max(positions)) return(NA_real_)
  hit <- match(coord, positions)
  if (!is.na(hit)) return(as.numeric(hit))
  approx(x=positions, y=seq_along(positions), xout=coord, ties="ordered")$y
}

all_pairwise_pass <- function(mat, idx, threshold) {
  z <- mat[idx, idx, drop=FALSE]
  v <- z[upper.tri(z)]
  length(v) > 0 && !any(is.na(v)) && all(v >= threshold)
}

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) return(data.frame(start_index=integer(), end_index=integer()))
  ord <- order(starts, ends)
  starts <- starts[ord]; ends <- ends[ord]
  cs <- starts[1]; ce <- ends[1]
  os <- integer(); oe <- integer()
  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= ce) {
        ce <- max(ce, ends[i])
      } else {
        os <- c(os, cs); oe <- c(oe, ce)
        cs <- starts[i]; ce <- ends[i]
      }
    }
  }
  os <- c(os, cs); oe <- c(oe, ce)
  data.frame(start_index=os, end_index=oe)
}

# Read r2 matrix
dat <- read.table(matrix_file, header=TRUE, sep="\t", check.names=FALSE,
                  stringsAsFactors=FALSE, quote="", comment.char="")
row_positions <- suppressWarnings(as.numeric(dat[[1]]))
if (any(is.na(row_positions))) stop("First column of matrix must be genomic positions.")
mat <- as.matrix(dat[, -1, drop=FALSE])
storage.mode(mat) <- "numeric"
col_positions <- suppressWarnings(as.numeric(colnames(mat)))
if (any(is.na(col_positions))) stop("Matrix column names must be genomic positions.")
if (nrow(mat) != ncol(mat)) stop("Matrix must be square.")
if (!setequal(row_positions, col_positions)) stop("Row and column SNP sets differ.")

ord <- order(row_positions)
positions <- row_positions[ord]
mat <- mat[ord, , drop=FALSE]
mat <- mat[, match(positions, col_positions), drop=FALSE]
rownames(mat) <- colnames(mat) <- as.character(positions)

# Read SNP list
snp <- read.table(snp_file, header=FALSE, stringsAsFactors=FALSE)
if (ncol(snp) < 2) stop("SNPlist.txt must have: position MAF")
snp <- snp[,1:2]
colnames(snp) <- c("Position","MAF")
snp$Position <- suppressWarnings(as.numeric(snp$Position))
snp$MAF <- suppressWarnings(as.numeric(snp$MAF))
if (any(is.na(snp$Position)) || any(is.na(snp$MAF))) stop("Non-numeric Position/MAF in SNPlist.")
if (anyDuplicated(snp$Position)) stop("Duplicate positions in SNPlist.")
missing <- setdiff(positions, snp$Position)
if (length(missing) > 0) stop("Missing SNP(s) in SNPlist: ", paste(missing, collapse=", "))
maf <- snp$MAF[match(positions, snp$Position)]
if (any(maf > 0.5)) warning("Some MAF values are >0.5; values are used exactly as supplied.")

# Anchors
anchors <- read_anchor(anchor_arg)
anchor_idx <- sapply(anchors, coord_to_index, positions=positions)
if (any(is.na(anchor_idx))) {
  warning("Anchor(s) outside plotted range ignored: ", paste(anchors[is.na(anchor_idx)], collapse=", "))
}
anchors_plot <- data.frame(coordinate=anchors[!is.na(anchor_idx)], index=anchor_idx[!is.na(anchor_idx)])

# Find qualifying 3-SNP windows
n <- length(positions)
triplets <- list()
k <- 0
if (n >= MIN_BLOCK_SNPS) {
  for (i in seq_len(n - MIN_BLOCK_SNPS + 1)) {
    idx <- i:(i + MIN_BLOCK_SNPS - 1)
    if (all_pairwise_pass(mat, idx, R2_THRESHOLD)) {
      k <- k + 1
      z <- mat[idx, idx, drop=FALSE]
      v <- z[upper.tri(z)]
      triplets[[k]] <- data.frame(
        triplet_id=k,
        start_index=min(idx),
        end_index=max(idx),
        start_position=positions[min(idx)],
        middle_position=positions[idx[2]],
        end_position=positions[max(idx)],
        min_pairwise_r2=min(v),
        mean_pairwise_r2=mean(v)
      )
    }
  }
}
triplets_df <- if (length(triplets)) do.call(rbind, triplets) else data.frame()

# Merge overlapping triplets and choose lowest-MAF marker
blocks_df <- data.frame()
if (nrow(triplets_df) > 0) {
  ints <- merge_intervals(triplets_df$start_index, triplets_df$end_index)
  blocks <- vector("list", nrow(ints))
  for (b in seq_len(nrow(ints))) {
    s <- ints$start_index[b]; e <- ints$end_index[b]; idx <- s:e
    m <- maf[idx]
    minmaf <- min(m)
    rep_idx <- idx[which(m == minmaf)[1]]
    z <- mat[idx, idx, drop=FALSE]
    v <- z[upper.tri(z)]
    v <- v[!is.na(v)]
    blocks[[b]] <- data.frame(
      block_id=paste0("LD", b),
      start_index=s,
      end_index=e,
      start_position=positions[s],
      end_position=positions[e],
      n_markers=length(idx),
      representative_index=rep_idx,
      representative_position=positions[rep_idx],
      representative_MAF=maf[rep_idx],
      min_observed_r2=ifelse(length(v), min(v), NA_real_),
      mean_observed_r2=ifelse(length(v), mean(v), NA_real_)
    )
  }
  blocks_df <- do.call(rbind, blocks)
}


# ----------------------------
# Inter-cluster r2 summary
# ----------------------------

intercluster_df <- data.frame()

if (nrow(blocks_df) >= 2) {

  intercluster_rows <- list()
  zrow <- 0

  for (a in 1:(nrow(blocks_df) - 1)) {
    for (b in (a + 1):nrow(blocks_df)) {

      idx_a <- blocks_df$start_index[a]:blocks_df$end_index[a]
      idx_b <- blocks_df$start_index[b]:blocks_df$end_index[b]

      cross_mat <- mat[idx_a, idx_b, drop=FALSE]

      n_possible <- length(cross_mat)
      observed_mask <- !is.na(cross_mat)
      n_observed <- sum(observed_mask)
      n_missing <- n_possible - n_observed

      if (n_observed > 0) {

        observed_values <- cross_mat[observed_mask]

        max_r2 <- max(observed_values)
        mean_r2 <- mean(observed_values)
        n_ge_threshold <- sum(observed_values >= R2_THRESHOLD)

        # First pair attaining the maximum r2.
        max_loc <- which(cross_mat == max_r2, arr.ind=TRUE)[1, ]

        max_pos_a <- positions[idx_a[max_loc[1]]]
        max_pos_b <- positions[idx_b[max_loc[2]]]

      } else {

        max_r2 <- NA_real_
        mean_r2 <- NA_real_
        n_ge_threshold <- 0
        max_pos_a <- NA_real_
        max_pos_b <- NA_real_
      }

      zrow <- zrow + 1

      intercluster_rows[[zrow]] <- data.frame(
        cluster1 = blocks_df$block_id[a],
        cluster2 = blocks_df$block_id[b],
        cluster1_start = blocks_df$start_position[a],
        cluster1_end = blocks_df$end_position[a],
        cluster1_n_markers = blocks_df$n_markers[a],
        cluster2_start = blocks_df$start_position[b],
        cluster2_end = blocks_df$end_position[b],
        cluster2_n_markers = blocks_df$n_markers[b],
        n_possible_pairs = n_possible,
        n_observed_pairs = n_observed,
        n_missing_pairs = n_missing,
        max_r2 = max_r2,
        max_r2_position_cluster1 = max_pos_a,
        max_r2_position_cluster2 = max_pos_b,
        mean_observed_r2 = mean_r2,
        n_pairs_ge_threshold = n_ge_threshold,
        linked_at_threshold = n_ge_threshold > 0
      )
    }
  }

  intercluster_df <- do.call(rbind, intercluster_rows)
}

# Long heatmap
heat <- expand.grid(y=seq_len(n), x=seq_len(n))
heat$r2 <- as.vector(mat)

axis_labels <- rep("", n)
lab_idx <- seq(1, n, by=LABEL_EVERY)
axis_labels[lab_idx] <- format(positions[lab_idx], scientific=FALSE, trim=TRUE)

p <- ggplot(heat, aes(x=x, y=y, fill=r2)) +
  geom_tile() +
  scale_fill_gradientn(
    colours=c("#FFFFFF","#FFF3C4","#F4A261","#D1495B","#7F0000"),
    values=scales::rescale(c(0,0.2,0.5,0.8,1)),
    limits=c(0,1),
    na.value="grey85",
    name=expression(r^2)
  ) +
  scale_x_continuous(breaks=seq_len(n), labels=axis_labels, expand=c(0,0)) +
  scale_y_reverse(breaks=seq_len(n), labels=axis_labels, expand=c(0,0)) +
  coord_fixed() +
  labs(
    x="Genomic position (GRCh38)",
    y="Genomic position (GRCh38)",
    title=expression(paste("Pairwise ", r^2, " heatmap with LD-block representatives"))
  ) +
  theme_minimal(base_size=10) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6),
    axis.text.y=element_text(size=6),
    plot.title=element_text(hjust=0.5)
  )

# Solid deletion/anchor line(s)
if (nrow(anchors_plot) > 0) {
  p <- p +
    geom_vline(data=anchors_plot, aes(xintercept=index), inherit.aes=FALSE,
               linewidth=0.65, linetype="solid", colour="black") +
    geom_hline(data=anchors_plot, aes(yintercept=index), inherit.aes=FALSE,
               linewidth=0.65, linetype="solid", colour="black")
}

# Dotted lowest-MAF representative line(s)
if (nrow(blocks_df) > 0) {
  p <- p +
    geom_vline(data=blocks_df, aes(xintercept=representative_index), inherit.aes=FALSE,
               linewidth=0.55, linetype="dotted", colour="black") +
    geom_hline(data=blocks_df, aes(yintercept=representative_index), inherit.aes=FALSE,
               linewidth=0.55, linetype="dotted", colour="black")
}

base <- sub("\\.[^.]+$", "", matrix_file)
png_file <- paste0(base, ".ldblock_heatmap.png")
pdf_file <- paste0(base, ".ldblock_heatmap.pdf")
block_file <- paste0(base, ".ldblocks.tsv")
triplet_file <- paste0(base, ".qualifying_triplets.tsv")
intercluster_file <- paste0(base, ".intercluster_r2.tsv")

ggsave(png_file, p, width=10.4, height=9.6, dpi=250)
ggsave(pdf_file, p, width=10.4, height=9.6)

write.table(blocks_df, block_file, sep="\t", row.names=FALSE, quote=FALSE)
write.table(triplets_df, triplet_file, sep="\t", row.names=FALSE, quote=FALSE)

write.table(
  intercluster_df,
  intercluster_file,
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

cat("R2 threshold:", R2_THRESHOLD, "\n")
cat("Minimum consecutive SNPs:", MIN_BLOCK_SNPS, "\n")
cat("Number of SNPs:", n, "\n")
cat("Qualifying 3-SNP windows:", nrow(triplets_df), "\n")
cat("Merged LD-block candidates:", nrow(blocks_df), "\n")
cat("Inter-cluster comparisons:", nrow(intercluster_df), "\n")
if (nrow(intercluster_df) > 0) {
  cat(
    "Cluster pairs with >=1 inter-cluster pair at r2 threshold:",
    sum(intercluster_df$linked_at_threshold),
    "\n"
  )
}
if (nrow(blocks_df) > 0) {
  print(blocks_df[,c("block_id","start_position","end_position","n_markers",
                     "representative_position","representative_MAF")], row.names=FALSE)
}
cat("Outputs:\n", png_file, "\n", pdf_file, "\n", block_file, "\n", triplet_file, "\n", intercluster_file, "\n")
