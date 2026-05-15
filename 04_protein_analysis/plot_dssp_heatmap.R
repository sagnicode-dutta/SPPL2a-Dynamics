library(ggplot2)
library(dplyr)

# Standard DSSP Color Scheme with Full Names
dssp_colors <- c(
  "H" = "#E11839", # Alpha Helix
  "G" = "#B14399", # 3-10 Helix
  "I" = "#711839", # Pi Helix
  "E" = "#F9F505", # Beta Sheet
  "B" = "#948E05", # Beta Bridge
  "T" = "#00AEEF", # Turn
  "S" = "#0072BC", # Bend
  "~" = "#FFFFFF"  # Coil
)

dssp_labels <- c(
  "H" = "Alpha Helix",
  "G" = "3-10 Helix",
  "I" = "Pi Helix",
  "E" = "Beta Sheet",
  "B" = "Beta Bridge",
  "T" = "Turn",
  "S" = "Bend",
  "~" = "Coil"
)

plot_heatmap <- function(csv_path, system_name, output_pdf) {
  df <- read.csv(csv_path)
  df$Structure <- factor(df$Structure, levels = names(dssp_colors))

  p <- ggplot(df, aes(x = Time, y = Residue, fill = Structure)) +
    geom_raster() +
    scale_fill_manual(values = dssp_colors, labels = dssp_labels, name = "Secondary Structure") +
    scale_y_continuous(breaks = seq(350, 400, 5)) +
    # Mark TM6a region - SOLID BLACK for bounds
    geom_hline(yintercept = c(368, 379), linetype = "solid", color = "black", linewidth = 1.5) +
    labs(title = paste("Secondary Structure Evolution:", system_name),
         subtitle = "TM6a Region (Residues 350-400) | Solid: TM6a (368-379)",
         x = "Time (ns)", y = "Residue Number") +
    theme_minimal(base_size = 30) + 
    theme(
      text = element_text(family = "sans", face = "bold"),
      plot.title = element_text(size = 34),
      plot.subtitle = element_text(size = 24, color = "black"),
      axis.text.x = element_text(size = 28, color = "black", face = "bold"),
      axis.text.y = element_text(size = 28, color = "black", face = "bold"),
      axis.title = element_text(size = 32, color = "black", face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 26)
    )
  
  ggsave(output_pdf, p, width = 20, height = 12)
}

plot_heatmap("analysis/9K92/dssp_tidy.csv", "Apo SPPL2a", "analysis/plots/tm6a_heatmap_9K92.pdf")
plot_heatmap("analysis/9K93/dssp_tidy.csv", "Bound SPPL2a", "analysis/plots/tm6a_heatmap_9K93.pdf")

print("DSSP Heatmaps updated with Black Dashed focus markers.")