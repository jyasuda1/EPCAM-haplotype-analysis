# plot_recombination_rate.R
# Input: map.txt
#   Column 1: physical position (bp)
#   Column 2: recombination rate (cM/Mb)
#
# Space- or tab-delimited files are accepted.
# The script automatically detects whether the first line is a header.

library(ggplot2)
library(scales)

input_file <- "W1510Strictpossiblemap.txt"
output_png <- "recombination_rate_plotSep4.png"
output_pdf <- "recombination_rate_plotsep4.pdf"

# Set TRUE only when an additional LOESS-smoothed curve is required.
add_loess <- FALSE
loess_span <- 0.05

# Detect whether the first two fields of the first line are numeric.
first_line <- readLines(input_file, n = 1, warn = FALSE)
first_fields <- strsplit(trimws(first_line), "\\s+")[[1]]

is_number <- function(x) {
  grepl("^[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?$", x)
}

has_header <- length(first_fields) < 2 ||
  !is_number(first_fields[1]) ||
  !is_number(first_fields[2])

map_data <- read.table(
  input_file,
  header = has_header,
  sep = "",
  stringsAsFactors = FALSE,
  comment.char = "#",
  check.names = FALSE
)

if (ncol(map_data) < 2) {
  stop("map.txt must contain at least two columns.")
}

# Use the first two columns regardless of their original column names.
map_data <- map_data[, 1:2]
colnames(map_data) <- c("position_bp", "recombination_rate")

map_data$position_bp <- suppressWarnings(as.numeric(map_data$position_bp))
map_data$recombination_rate <- suppressWarnings(
  as.numeric(map_data$recombination_rate)
)

map_data <- map_data[
  is.finite(map_data$position_bp) &
  is.finite(map_data$recombination_rate),
]

if (nrow(map_data) == 0) {
  stop("No valid numeric data were found in map.txt.")
}

# Average duplicated physical positions, if present, and sort by position.
map_data <- aggregate(
  recombination_rate ~ position_bp,
  data = map_data,
  FUN = mean
)
map_data <- map_data[order(map_data$position_bp), ]

# Convert the genomic coordinate from bp to Mb for plotting.
map_data$position_Mb <- map_data$position_bp / 1e6

p <- ggplot(
  map_data,
  aes(x = position_Mb, y = recombination_rate)
) +
  geom_line(linewidth = 0.55) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.001)
  ) +
  labs(
    title = "Recombination rate around the core haplotype region",
    x = "Physical position (Mb)",
    y = "Recombination rate (cM/Mb)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Optional LOESS curve. The raw line remains visible underneath it.
if (add_loess) {
  p <- p +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = FALSE,
      linewidth = 0.8
    )
}

ggsave(
  filename = output_png,
  plot = p,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 9,
  height = 5.5
)

print(p)

cat("Input points:", nrow(map_data), "\n")
cat(
  "Physical range:",
  sprintf("%.6f–%.6f Mb", min(map_data$position_Mb), max(map_data$position_Mb)),
  "\n"
)
cat("Saved:", output_png, "and", output_pdf, "\n")
