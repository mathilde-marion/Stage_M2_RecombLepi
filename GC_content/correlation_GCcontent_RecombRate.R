# ═════════════════════════════════════════════════════════════════
# Corrélation recomb relatif ~ GC content pour toutes les espèces ####
# ═════════════════════════════════════════════════════════════════

# A FAIRE TOURNER AVANT : ~/Bureau/DataPhenolep/recombRateRelative_GCcontent/GC_recombRateRelativeRelative/script_correlation_GCcontent_recombRateRelativeRelative.R
library(tidyverse)
library(ggplot2)
library(patchwork)

species_colors <- c(
  "lys_bel.POP1" = "#397df4", "lys_bel.POP2" = "#7ea8f6",
  "mel_pho.POP1" = "#38a6a5", "mel_pho.POP2" = "#7fc8c7",
  "pie_nap.POP1" = "#fff52e", "pie_nap.POP2" = "#fff87a",
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

# ── Chargement des données ────────────────────────────────────────────────────

data_dir <- "~/Bureau/DataPhenolep/RecombRate_GCcontent/GC_recombRateRelative"

files <- list.files(
  data_dir,
  pattern = "\\.RecombRateRelative_GCcontent\\.txt$",
  full.names = TRUE
)

recomb_gc_all <- bind_rows(lapply(files, function(f) {
  df <- read.table(f, header = TRUE)
  df$chrom <- as.character(df$chrom)
  base <- gsub("\\.RecombRateRelative_GCcontent\\.txt$", "", basename(f))
  df$species <- gsub("\\..*", "", base)
  df$pop     <- gsub(".*\\.", "", base)
  df
}))

# ── Nettoyage ─────────────────────────────────────────────────────────────────

recomb_gc_clean <- recomb_gc_all %>%
  filter(
    recombRateRelative > 0,
    is.finite(recombRateRelative),
    !is.na(GC)
  ) %>%
  mutate(
    GC_percent = GC * 100,
    sp_pop     = paste(species, pop, sep = ".")  # clé pour la palette
  )


# ── Graphiques par espèce / population ───────────────────────────────────────

plot_list <- list()

for (sp in unique(recomb_gc_clean$species)) {
  for (pop in unique(recomb_gc_clean$pop[recomb_gc_clean$species == sp])) {
    
    d <- recomb_gc_clean %>%
      filter(species == sp, pop == pop)
    
    sp_pop_key <- paste(sp, pop, sep = ".")
    line_color <- species_colors[sp_pop_key]
    # fallback si l'espèce n'est pas dans la palette
    if (is.na(line_color)) line_color <- "black"
    
    ct <- cor.test(d$recombRateRelative, d$GC_percent, method="spearman")
    
    r <- round(ct$estimate, 3)
    p <- ct$p.value
    
    label_r <- ifelse(
      p < 0.001,
      paste0("r = ", r, "\np < 0.001"),
      paste0("r = ", r, "\np = ", signif(p, 3))
    )
    
    plot_list[[sp_pop_key]] <-
      ggplot(d, aes(x = recombRateRelative, y = GC_percent)) +
      geom_point(alpha = 0.15, size = 0.8, color = "grey20") +
      geom_smooth(method = "lm", color = line_color, se = FALSE, linewidth = 0.9) +
      annotate("text", x = -Inf, y = Inf, label = label_r,
               hjust = -0.1, vjust = 1.1, size = 3.5) +
      labs(
#        title = sp_pop_key,
        x = "Taux de recombinaison relatif",
        y = "Contenu en GC (%)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        text              = element_text(family = "Montserrat", color = "grey20"),
        axis.title        = element_text(face = "bold", size = 12, color = "grey20"),
        axis.text         = element_text(size = 10,  color = "grey20"),
        axis.line         = element_line(color = "grey20"),
        axis.ticks        = element_line(color = "grey20"),
        legend.text       = element_text(size = 10, color = "grey20"),
        legend.position   = "right",
        legend.key        = element_rect(fill = NA, color = NA),
        panel.background  = element_rect(fill = "white", color = NA),
        plot.background   = element_rect(fill = "white", color = NA),
        panel.grid        = element_blank()
      )
  }
}

# ── Assemblage ───────────────────────────────────────────────────────────────

final_plot <- wrap_plots(plot_list, ncol = 4)

print(final_plot)

ggsave(
  file.path(data_dir, "recombRateRelative_GC_by_pop_presque_final_pr_annexe.pdf"),
  final_plot,
  width  = 4.5 * ceiling(sqrt(length(plot_list))),
  height = 4 * ceiling(sqrt(length(plot_list))),
  device = cairo_pdf)

# ═════════════════════════════════════════════════════════════════
# Corrélation taux recomb ~ GC content pour toutes les espèces ####
# ═════════════════════════════════════════════════════════════════

# A FAIRE TOURNER AVANT : :~/Bureau/DataPhenolep/recombRateRelative_GCcontent/GC_recombRateRelative/script_correlation_GCcontent_recombRateRelative.R
library(tidyverse)
library(ggplot2)
library(patchwork)

species_colors <- c(
  "lys_bel.POP1" = "#397df4", "lys_bel.POP2" = "#7ea8f6",
  "mel_pho.POP1" = "#38a6a5", "mel_pho.POP2" = "#7fc8c7",
  "pie_nap.POP1" = "#fff52e", "pie_nap.POP2" = "#fff87a",
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


# ── Chargement ────────────────────────────────────────────────────────────────

data_dir <- "~/Bureau/DataPhenolep/RecombRate_GCcontent/GC_recombRate"

files <- list.files(data_dir, pattern = "*.RecombRate_GCcontent.txt", full.names = TRUE)

recomb_gc_all <- do.call(rbind, lapply(files, function(f) {
  df <- read.table(f, header = TRUE)
  base <- gsub("\\.RecombRate_GCcontent\\.txt$", "", basename(f))  # ex: "lys_bel.POP1"
  df$species <- gsub("\\..*", "", base)   # "lys_bel"
  df$pop     <- gsub(".*\\.", "", base)   # "POP1"
  df
}))

# ── Nettoyage ─────────────────────────────────────────────────────────────────

recomb_gc_clean <- recomb_gc_all[
  recomb_gc_all$recombRate > 0 &
    is.finite(recomb_gc_all$recombRate) &
    !is.na(recomb_gc_all$GC), ]


recomb_gc_clean <- recomb_gc_clean %>%
  mutate(
    recombRate_cMMb  = recombRate * 1e8,
    GC_percent = GC * 100,
    sp_pop     = paste(species, pop, sep = ".")
  )

# ── Graphique par population ──────────────────────────────────────────────────

plot_list <- list()

for (sp in unique(recomb_gc_clean$species)) {
  for (pop in unique(recomb_gc_clean$pop[recomb_gc_clean$species == sp])) {
    
    d <- recomb_gc_clean %>%
      filter(species == sp, pop == pop)
    
    sp_pop_key <- paste(sp, pop, sep = ".")
    line_color <- species_colors[sp_pop_key]
    if (is.na(line_color)) line_color <- "grey30"
    
    ct <- cor.test(d$recombRate_cMMb, d$GC_percent, method = "spearman")
    
    r <- round(ct$estimate, 3)
    p <- ct$p.value
    
    label_r <- ifelse(
      p < 0.001,
      paste0("r = ", r, "\np < 0.001"),
      paste0("r = ", r, "\np = ", signif(p, 3))
    )
    
    plot_list[[sp_pop_key]] <-
      ggplot(d, aes(x = recombRate_cMMb, y = GC_percent)) +
      
      geom_point(alpha = 0.15, size = 0.8, color = "grey20") +
      
      geom_smooth(
        method = "lm",
        color = line_color,
        se = FALSE,
        linewidth = 0.9
      ) +
      
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = label_r,
        hjust = -0.1, vjust = 1.1,
        size = 3.5
      ) +
      
      labs(
        y = "Contenu en GC (%)",
        x = "Taux de recombinaison (cM/Mb)"
      ) +
      
      theme_classic(base_size = 12) +
      
      theme(
        text = element_text(family = "Montserrat", color = "grey20"),
        
        axis.title = element_text(face = "bold", size = 12, color = "grey20"),
        axis.text  = element_text(size = 10, color = "grey20"),
        
        axis.line  = element_line(color = "grey20"),
        axis.ticks = element_line(color = "grey20"),
        
        plot.title = element_blank(),
        legend.position = "none",
        
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
      )
  }
}
# ── Assemblage patchwork ──────────────────────────────────────────────────────

final_plot <- wrap_plots(plot_list, ncol=4)

print(final_plot)

ggsave(
  file.path(data_dir, "recombRate_vs_GC_by_pop_presque_final_pr_annexe_ptet.pdf"),
  final_plot,
  width  = 4.5 * ceiling(sqrt(length(plot_list))),
  height = 4 * ceiling(sqrt(length(plot_list))),
  device = cairo_pdf)

# ════════════════════════════════════════════════════
# Corrélation LOG taux recomb ~ GC content pour toutes les espèces
# ════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Chargement ────────────────────────────────────────────────────────────────

data_dir <- "~/Bureau/DataPhenolep/RecombRate_GCcontent/GC_recombRate"

files <- list.files(data_dir, pattern = "*.RecombRate_GCcontent.txt", full.names = TRUE)

recomb_gc_all <- do.call(rbind, lapply(files, function(f) {
  df <- read.table(f, header = TRUE)
  base <- gsub("\\.RecombRate_GCcontent\\.txt$", "", basename(f))
  df$species <- gsub("\\..*", "", base)
  df$pop     <- gsub(".*\\.", "", base)
  df
}))

# ── EXCLUSION ──
recomb_gc_all <- recomb_gc_all %>%
  mutate(sp_pop = paste(species, pop, sep = ".")) %>%
  filter(
    !sp_pop %in% c(
      "pap_hos.POP1",
      "pap_hos.POP2",
      "mel_pho.POP1",
      "mel_gal.POP1"
    )
  )

# ── Nettoyage ─────────────────────────────────────────────────────────────────

recomb_gc_clean <- recomb_gc_all[
  recomb_gc_all$recombRate > 0 &
    is.finite(recomb_gc_all$recombRate) &
    !is.na(recomb_gc_all$GC), ]

recomb_gc_clean <- recomb_gc_clean %>%
  mutate(
    log_recomb = log10(recombRate),
    GC_percent = GC * 100,
    sp_pop     = paste(species, pop, sep = ".")
  )

# ── Graphiques ────────────────────────────────────────────────────────────────
plot_list <- list()

for (sp_pop_key in unique(recomb_gc_clean$sp_pop)) {
  
  d <- recomb_gc_clean %>%
    filter(sp_pop == sp_pop_key)
  
  line_color <- species_colors[sp_pop_key]
  if (is.na(line_color)) line_color <- "grey30"
  
  ct <- cor.test(d$log_recomb, d$GC_percent, method = "spearman")
  r  <- round(ct$estimate, 3)
  p  <- ct$p.value
  
  label_r <- ifelse(
    p < 0.001,
    paste0("r = ", r, "\np < 0.001"),
    paste0("r = ", r, "\np = ", signif(p, 3))
  )
  
  plot_list[[sp_pop_key]] <-
    ggplot(d, aes(x = log_recomb, y = GC_percent)) +
    
  #  geom_point(alpha = 0.15, size = 0.8, color = "grey20") +
    geom_point_rast(
      alpha = 0.15,
      size = 0.8,
      color = "grey20",
      raster.dpi = 300
    ) +
    
    geom_smooth(
      method = "lm",
      color = line_color,
      se = FALSE,
      linewidth = 0.9
    ) +
    
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_r,
      hjust = -0.1, vjust = 1.1,
      size = 7
    ) +
    
    labs(
#      title = sp_pop_key,
      y = "Contenu en GC (%)",
      x = "Logarithme du taux de recombinaison"
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      text = element_text(family = "Montserrat", color = "grey20"),
      axis.title = element_text(face = "bold", size = 19),
      axis.text  = element_text(size = 17),
      legend.position = "none"
    )
}

# ── Assemblage ────────────────────────────────────────────────────────────────

final_plot <- wrap_plots(plot_list, ncol = 4)

print(final_plot)

ggsave(
  file.path("~/Bureau/latex/figures/annexes/GC_recombLog10_by_pop_annexe_rastr_test.pdf"),
  final_plot,
  width  = 33,
  height = 22,
  device = cairo_pdf)
