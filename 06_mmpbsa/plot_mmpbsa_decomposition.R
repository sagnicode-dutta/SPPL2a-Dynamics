library(ggplot2)
library(dplyr)
library(tidyr)

# Load SPPL2a data
sppl2a_data <- read.table("analysis/mmpbsa/sppl2a_mmpbsa/top_residues_sppl2a_final.txt", 
                          col.names = c("Residue", "Energy", "SEM"))
sppl2a_data$System <- "Bound SPPL2a"

# Load PS1 data
ps1_data <- read.table("analysis/mmpbsa/ps1_mmpbsa/top_residues_ps1_final.txt", 
                       col.names = c("Residue", "Energy", "SEM"))
ps1_data$System <- "PS1 (Bound)"

# Combine data
combined_data <- rbind(sppl2a_data, ps1_data)

# Extract residue name and number
combined_data <- combined_data %>%
  mutate(ResName = sub(".*:.*:(.*):.*", "\\1", Residue),
         ResNum = sub(".*:.*:.*:(.*)", "\\1", Residue)) %>%
  mutate(ResLabel = paste0(ResName, ResNum))

# APPLY EXCLUSIONS AND INCLUSIONS
# 1. PS1 Exclusions: ARG269, GLU273
# 2. Threshold: Energy < -20 kcal/mol
# 3. Always include Catalytic Dyads even if they were slightly different (though here they are strong)

plot_data <- combined_data %>%
  # Explicitly exclude artifacts
  filter(!(System == "PS1 (Bound)" & ResLabel %in% c("ARG269", "GLU273"))) %>%
  # Filter by energy threshold
  filter(Energy < -20) %>%
  # Exclude ligand itself if needed (FTO 495 in SPPL2a)
  filter(!grepl("FTO", ResLabel))

# Define colors
system_colors <- c("Bound SPPL2a" = "#53c719", "PS1 (Bound)" = "#9929ea")

# Create a combined label for the Y axis to distinguish systems if residues have same number
plot_data <- plot_data %>%
  mutate(PlotLabel = paste0(ResLabel, " (", ifelse(System == "Bound SPPL2a", "SPPL2a", "PS1"), ")"))

# Create plot
p <- ggplot(plot_data, aes(x = reorder(PlotLabel, -Energy), y = Energy, fill = System)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = Energy - SEM, ymax = Energy + SEM), width = 0.2) +
  coord_flip() +
  scale_fill_manual(values = system_colors) +
  theme_bw() +
  labs(title = "Top Binding Energy Contributors (MMPBSA)",
       subtitle = "Excluded artifacts (ARG269, GLU273); Threshold < -20 kcal/mol",
       x = "Residue",
       y = "Decomposition Energy (kcal/mol)",
       fill = "System") +
  theme(axis.text.y = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

# Save plot
ggsave("analysis/plots/mmpbsa_decomposition_refined.pdf", p, width = 8, height = 7)
ggsave("analysis/plots/mmpbsa_decomposition_refined.png", p, width = 8, height = 7)
