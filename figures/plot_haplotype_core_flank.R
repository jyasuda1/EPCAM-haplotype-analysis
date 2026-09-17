# Haplotype core/flank plot
# Core region is calculated as the intersection of all case intervals:
# max(start) to min(end)

library(ggplot2)
library(dplyr)
library(tidyr)

dat <- read.table(text = "Cases	start_cM	end_cM	ColorCore	Colorflank


", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA")

deletion_cM <- 69.1352076322698

core_start <- max(dat$start_cM)
core_end   <- min(dat$end_cM)

# Preserve the input order from top to bottom
dat$Cases <- factor(dat$Cases, levels = rev(dat$Cases))

segments <- dat %>%
  rowwise() %>%
  do(data.frame(
    Cases = .$Cases,
    part = c("left", "core", "right"),
    x = c(.$start_cM, core_start, core_end),
    xend = c(core_start, core_end, .$end_cM),
    color_group = c(.$Colorflank, "core", .$Colorflank)
  )) %>%
  ungroup() %>%
  filter(xend > x) %>%
  mutate(
    color_group = ifelse(is.na(color_group), "core", color_group)
  )

color_values <- c(
  core = "black",
  light_green = "lightgreen",
  orange = "orange",
  red = "red"
)

p <- ggplot(segments) +
  geom_segment(
    aes(x = x, xend = xend, y = Cases, yend = Cases, color = color_group),
    linewidth = 3.5,
    lineend = "butt"
  ) +
  geom_vline(
    xintercept = deletion_cM,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  scale_color_manual(values = color_values, guide = "none") +
  labs(
    title = "Shared haplotype intervals",
    subtitle = sprintf("Core: %.6f–%.6f cM; deletion: %.6f cM",
                       core_start, core_end, deletion_cM),
    x = "Genetic position (cM)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "haplotype_core_flank_plotSept4.png",
  plot = p,
  width = 10,
  height = 6.5,
  dpi = 300
)

print(p)
