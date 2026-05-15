library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Load data
apo_df <- read.csv("analysis/dccm/dccm_apo.csv")
bound_df <- read.csv("analysis/dccm/dccm_bound.csv")

# 1. Identify present residues and map to 1..N indices
res_present <- sort(unique(c(apo_df$Residue_i, bound_df$Residue_i)))
n_res <- length(res_present)
res_map <- setNames(1:n_res, res_present)

# 2. TM Definitions
tm_defs <- data.frame(
  name = c("TM2", "TM6", "TM6a", "TM9"),
  start = c(1, 22, 46, 58),
  end = c(21, 45, 57, 75),
  color = c("#FF8080", "#80FFE5", "#80D4FF", "#A9FFD1")
)

# 3. Component Factory Functions
create_heatmap <- function(df) {
  df$idx_i <- res_map[as.character(df$Residue_i)]
  df$idx_j <- res_map[as.character(df$Residue_j)]
  
  lbl_idx <- seq(1, n_res, by=5)
  lbl_val <- res_present[lbl_idx]
  
  ggplot(df, aes(x = idx_i, y = idx_j, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1), guide = "none") +
    scale_x_continuous(breaks = lbl_idx, labels = lbl_val, expand = c(0,0)) +
    scale_y_continuous(breaks = lbl_idx, labels = lbl_val, expand = c(0,0)) +
    coord_fixed() +
    labs(x = "Residue Number", y = "Residue Number") +
    theme_bw(base_size = 36) +
    theme(
      text = element_text(family = "sans", face = "bold", color = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 26, color = "black", face = "bold"),
      axis.text.y = element_text(size = 26, color = "black", face = "bold"),
      axis.title = element_text(size = 38, color = "black", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
      panel.grid = element_blank(),
      plot.margin = margin(t=0, r=10, b=10, l=10) # No top margin to stay close to title
    )
}

create_left_bar <- function() {
  ggplot(tm_defs) +
    geom_rect(aes(ymin = start - 0.5, ymax = end + 0.5, xmin = 0, xmax = 1, fill = color), 
              color = "black", linewidth = 1.5) +
    geom_text(aes(y = (start + end)/2, x = 0.5, label = name), 
              family = "sans", fontface = "bold", color = "black", size = 12, angle = 90) +
    scale_fill_identity() +
    # Lock limits and expand to match heatmap Y-axis exactly
    scale_y_continuous(limits = c(0.5, n_res + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t=0, r=5, b=10, l=20)) # Added left padding (l=20)
}

# 4. Assembly of a Single Panel
assemble_system_panel <- function(df, title_str) {
  p_main <- create_heatmap(df)
  p_left <- create_left_bar()
  
  # Heading object with very tight bottom margin to move it closer to plot
  p_title <- ggplot() + 
    annotate("text", x = 0, y = 0, label = title_str, family = "sans", fontface = "bold", color = "black", size = 22) + 
    theme_void() + theme(plot.margin = margin(t=10, b=-20))

  # Use align and axis to force bars to match plot area
  lower_row <- plot_grid(p_left, p_main, ncol = 2, rel_widths = c(0.07, 1), align = "h", axis = "bt")
  
  # Center title only over the main plot box
  upper_row_title <- plot_grid(NULL, p_title, ncol = 2, rel_widths = c(0.07, 1))
  
  # Combine Title and Plot
  final <- plot_grid(upper_row_title, lower_row, ncol = 1, rel_heights = c(0.05, 1))
  return(final)
}

# 5. Generate Panels
panel_apo <- assemble_system_panel(apo_df, "Apo SPPL2a")
panel_bound <- assemble_system_panel(bound_df, "Bound SPPL2a")

# Shared Legend with extra title spacing
p_leg_dummy <- ggplot(apo_df, aes(x=Residue_i, y=Residue_j, fill=Correlation)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-1,1), name="Correlation") +
  theme(legend.title=element_text(family="sans", face="bold", color="black", size=32, margin=margin(b=30)),
        legend.text=element_text(family="sans", face="bold", color="black", size=28),
        legend.key.height=unit(5, "cm"))
leg <- get_legend(p_leg_dummy)

# Final Grid Assembly
# Add a leading NULL for overall left padding
final_output <- plot_grid(NULL, panel_apo, NULL, panel_bound, leg, 
                          ncol = 5, rel_widths = c(0.05, 1, 0.2, 1, 0.25))

# Save
ggsave("analysis/plots/tm2_tm6a_dccm_comparison.pdf", final_output, width = 42, height = 20, device = cairo_pdf)

print("DCCM comparison plot successfully regenerated with aligned rectangles, padded legend, and lowered titles.")
