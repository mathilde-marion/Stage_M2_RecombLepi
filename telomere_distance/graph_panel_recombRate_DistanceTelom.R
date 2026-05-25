library(tidyverse)
library(patchwork)

#--------------------------------------------------------------
# Palette de couleurs par espèce/population
#--------------------------------------------------------------

species_colors <- c(
  "lys_bel.POP1" = "#397df4", "lys_bel.POP2" = "#7ea8f6",
  "mel_pho.POP1" = "#38a6a5", "mel_pho.POP2" = "#7fc8c7",
  "pie_nap.POP1" = "#d4c000", "pie_nap.POP2" = "#e6d94a",
  "pap_mac.POP1" = "#ff5a80", "pap_mac.POP2" = "#ff9ab0",
  "pap_hos.POP1" = "#ffa7dd", "pap_hos.POP2" = "#ffd1ee",
  "ito_sal.POP1" = "#73af48", "ito_sal.POP2" = "#a8d38a",
  "hel_mel.POP1" = "#f73131", "hel_mel.POP2" = "#fa7b7b",
  "hel_cyd.POP1" = "#1d6996", "hel_cyd.POP2" = "#5c97bd",
  "hel_tim.POP1" = "#0f8554", "hel_tim.POP2" = "#5fbf9a",
  "hel_num.POP1" = "#edad08", "hel_num.POP2" = "#f3cc5c",
  "pol_ica.POP1" = "#8e38a4", "pol_ica.POP2" = "#b97bca",
  "ant_car.POP1" = "#784421", "ant_car.POP2" = "#a97a5a",
  "mel_gal.POP1" = "#ff672e", "mel_gal.POP2" = "#ff9a73"
)

#--------------------------------------------------------------
# Chargement de tous les fichiers
#--------------------------------------------------------------

data_dir <- "~/Bureau/DataPhenolep/OutputRelernnCluster/sliding_windows"
files <- list.files(
  data_dir,
  pattern = "\\.quantiles\\.SlidingWindow100k-1M\\.txt$",
  full.names = TRUE
)

all_data <- bind_rows(lapply(files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  df <- df %>% mutate(Chr = str_remove_all(Chr, "b'|'"))
  base       <- gsub("\\.quantiles\\.SlidingWindow100k-1M\\.txt$", "", basename(f))
  df$species <- sub("\\.POP[12]$", "", base)
  df$pop     <- sub(".*\\.", "", base)
  df$sp_pop  <- paste(df$species, df$pop, sep = ".")
  df
}))

all_data <- all_data %>%
  filter(
    !sp_pop %in% c(
      "pap_hos.POP1",
      "pap_hos.POP2",
      "mel_pho.POP1",
      "mel_gal.POP1"
    )
  )
cat("Espèces/populations détectées :\n")
print(unique(all_data$sp_pop))

#--------------------------------------------------------------
# Pipeline : fold + mean pour toutes les espèces/pop
#--------------------------------------------------------------

all_means <- bind_rows(lapply(unique(all_data$sp_pop), function(sp_pop_key) {
  
  data <- all_data %>% filter(sp_pop == sp_pop_key)
  
  chrom_sizes <- data %>%
    group_by(Chr) %>%
    summarise(chrom_length = max(Window_End), .groups = "drop")
  
  data <- data %>%
    left_join(chrom_sizes, by = "Chr") %>%
    mutate(
      mid_window     = (Window_Start + Window_End) / 2,
      dist_to_end    = pmin(mid_window, chrom_length - mid_window),
      dist_to_end_Mb = dist_to_end / 1e6
    )
  
  data_folded <- data %>%
    group_by(Chr, dist_to_end_Mb) %>%
    summarise(mean_quantile = mean(Mean_recomb_quantile_100, na.rm = TRUE),
              .groups = "drop")
  
  data_mean <- data_folded %>%
    group_by(dist_to_end_Mb) %>%
    summarise(
      mean_quantile = mean(mean_quantile, na.rm = TRUE),
      n_chr         = n(),
      .groups = "drop"
    ) %>%
    filter(n_chr >= 5) %>%
    mutate(sp_pop = sp_pop_key)
  
  data_mean
}))

#--------------------------------------------------------------
# Limites fixes des axes (identiques pour tous les panels)
#--------------------------------------------------------------

x_max <- max(all_means$dist_to_end_Mb, na.rm = TRUE)
y_lim <- c(30, 60)
x_lim <- c(0, x_max)

#--------------------------------------------------------------
# Graphique unique avec toutes les espèces/pop + légende
#--------------------------------------------------------------

plot_all <-
  ggplot(all_means, aes(x = dist_to_end_Mb, y = mean_quantile,
                        color = sp_pop, group = sp_pop)) +
  geom_line(linewidth = 0.8, na.rm = TRUE, alpha = 0.6) +
  scale_color_manual(values = species_colors) +
  scale_x_continuous(limits = x_lim, breaks = c(0, 2.5, 5, 7.5, 10, 12.5)) +
  scale_y_continuous(limits = y_lim, breaks = c(30, 40, 50, 60 )) +
  labs(
    x     = "Distance au télomère (Mb)",
    y     = "Taux de recombinaison relatif",
    color = NULL  
  ) +
  guides(
    color = guide_legend(
      keywidth  = unit(1.5, "cm"), 
      keyheight = unit(0.4, "cm")
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = "none",  
    text              = element_text(family = "Montserrat", color = "grey20"),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title        = element_text(face = "bold", size = 14, color = "grey20"),
    axis.text         = element_text(size = 13,  color = "grey20"),
    axis.line         = element_line(color = "grey20"),
    axis.ticks        = element_line(color = "grey20"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.grid        = element_blank()
  )

print(plot_all)

##
make_plot <- function(data, sp_pop_key) {
  
  line_color <- species_colors[sp_pop_key]
  if (is.na(line_color)) line_color <- "steelblue"
  
  # Taille chromosomes
  chrom_sizes <- data %>%
    group_by(Chr) %>%
    summarise(chrom_length = max(Window_End), .groups = "drop")
  
  # Distance au télomère
  data <- data %>%
    left_join(chrom_sizes, by = "Chr") %>%
    mutate(
      mid_window     = (Window_Start + Window_End) / 2,
      dist_to_end    = pmin(mid_window, chrom_length - mid_window),
      dist_to_end_Mb = dist_to_end / 1e6
    )
  
  # -----------------------------
  # NON MOYENNÉ (chromosomes)
  # -----------------------------
  data_folded <- data %>%
    group_by(Chr, dist_to_end_Mb) %>%
    summarise(
      mean_quantile = mean(Mean_recomb_quantile_100, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # MOYENNE
  # -----------------------------
  data_mean <- data_folded %>%
    group_by(dist_to_end_Mb) %>%
    summarise(
      mean_quantile = mean(mean_quantile, na.rm = TRUE),
      n_chr         = n(),
      .groups = "drop"
    ) %>%
    filter(n_chr >= 5)
  
  # -----------------------------
  # PLOT
  # -----------------------------
  ggplot() +
    
    # lignes chromosomes (fines)
    geom_line(
      data = data_folded,
      aes(x = dist_to_end_Mb, y = mean_quantile, group = Chr),
      color = line_color,
      linewidth = 0.3,
      alpha = 0.3
    ) +
    
    # ligne moyenne (épaisse)
    geom_line(
      data = data_mean,
      aes(x = dist_to_end_Mb, y = mean_quantile),
      color = line_color,
      linewidth = 1
    ) +
    
    labs(
      x = "Distance au télomère (Mb)",
      y = "Taux de recombinaison relatif"
    ) +
    
    scale_x_continuous(breaks = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)) +
    
    theme_classic(base_size = 12) +
    theme(
      axis.title        = element_text(face = "bold", size = 20, color = "grey20"),
      axis.text         = element_text(size = 18, color = "grey20"),
      axis.line         = element_line(color = "grey20"),
      axis.ticks        = element_line(color = "grey20"),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.grid        = element_blank()
    )
}

#--------------------------------------------------------------
# Génération de tous les graphiques
#--------------------------------------------------------------

plot_list <- list()

for (sp_pop_key in sort(unique(all_data$sp_pop))) {
  d <- all_data %>% filter(sp_pop == sp_pop_key)
  plot_list[[sp_pop_key]] <- make_plot(d, sp_pop_key)
}

#--------------------------------------------------------------
# Assemblage avec patchwork
#--------------------------------------------------------------

final_plot <- wrap_plots(plot_list, ncol = 4)
print(final_plot)

ggsave(
  file.path("~/Bureau/latex/figures/annexes/telomere_annexe.pdf"),
  final_plot,
  width  = 35,
  height = 24.74,
  device = cairo_pdf)


#pannel pour rapport
p_lys <- make_plot(
  all_data %>% filter(sp_pop == "lys_bel.POP2"),
  "lys_bel.POP2"
)

p_ito <- make_plot(
  all_data %>% filter(sp_pop == "hel_mel.POP1"),
  "hel_mel.POP1"
)

left_panel <- p_lys / p_ito
combined_plot <- left_panel | plot_all

ggsave(
  file.path("~/Bureau/latex/figures/telomere_recombRateRelative_VF.pdf"),
  combined_plot,
  width = 12,
  height = 6,
  device = cairo_pdf
)

