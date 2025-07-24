plot_strand_scatter <- function(ifn_tab, ifn_codon, title, fdir) {
  # Read the strand comparison table
  df <- read.delim(ifn_tab)

  # Read codon info
  ref <- read.delim(ifn_codon, colClasses = c(NULL, "character", NULL))
  ref_id <- paste(ref$coord, ref$aa, ref$codon, sep = "_")
  df$is_ref <- df$id == ref_id

  library(ggplot2)
  library(ggrepel)
  dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

  # Create labels by dropping string before first underscore
  df$label <- sub("^[^_]*_", "", df$id)

  # Create scatter plot with codon labels
  p <- ggplot(df, aes(x = plus, y = minus)) +
    geom_point(aes(color = is_ref), size = 3) +
    scale_color_manual(values = c("FALSE" = "#2c7bb6", "TRUE" = "#d7191c")) +
    geom_text_repel(aes(label = label),
      size = 3.5,
      box.padding = 0.5,
      point.padding = 0.2,
      force = 10
    ) +
    theme_minimal() +
    labs(
      title = sprintf("Position %s\nReference: %s", title, ref_id),
      x = "Plus Strand Log-Likelihood",
      y = "Minus Strand Log-Likelihood"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    coord_fixed() # Make plot square with equal aspect ratio

  # Save plot as PDF
  ofn <- file.path(fdir, sprintf("PvM_%s.pdf", title))
  ggsave(ofn, p, width = 8, height = 8)

  cat("Plot saved to:", ofn, "\n")
}
