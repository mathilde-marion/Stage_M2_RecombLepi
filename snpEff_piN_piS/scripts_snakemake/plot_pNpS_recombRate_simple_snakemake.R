#!/usr/bin/env Rscript
# Usage: Rscript plot_pNpS_recomb.R <input.tsv> <output.png>

library(ggplot2)
library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_pNpS_recomb.R <input.tsv> <output.png>")

in_tsv  <- args[1]
out_png <- args[2]

# chargement et nettoyage des données
df <- read.table(in_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df[, c("pN_pS", "recombRate")]
df <- df[df$pN_pS != "NA", ]
df$pN_pS      <- as.numeric(df$pN_pS)
df$recombRate <- as.numeric(df$recombRate)
df <- na.omit(df)

df <- df %>% mutate(recombRate_cMMb  = recombRate  * 1e8)

# stats
ct  <- cor.test(df$recombRate_cMMb, df$pN_pS, method = "spearman")
rho <- round(ct$estimate, 3)
pv  <- signif(ct$p.value, 2)
lab <- paste0("Spearman ρ = ", rho, "\np = ", pv)

# plot
p  <- ggplot(df, aes(x = recombRate_cMMb, y = pN_pS)) +
  geom_point(color = "#2c2c2c", size = 1.8, alpha = 0.55, shape = 16) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#c0392b", fill = "#c0392b", linewidth = 0.8, alpha = 0.15) +
  annotate("text",
           x = Inf, y = Inf,
           label = lab,
           hjust = 1.1, vjust = 1.5,
           size = 3.8, family = "mono", color = "#2c2c2c") +
  labs(
    x = "Taux de recombinaison (cM/Mb)",
    y = expression(bold(p[N]/p[S]))
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x  = element_text(face = "bold", size = 13),
    axis.title.y  = element_text(face = "bold", size = 13),
    axis.text     = element_text(color = "#2c2c2c"),
    axis.line     = element_line(color = "#2c2c2c", linewidth = 0.5),
    axis.ticks    = element_line(color = "#2c2c2c"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(out_png, plot = p, width = 6, height = 5, dpi = 300, bg = "white")
cat("Done ->", out_png, "\n")
