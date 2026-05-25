# ════════════════════════════════════════════════════════
# Calcul des quantiles (1–100) pour tous les fichiers ####
# ════════════════════════════════════════════════════════

# ── 1. Dossier ──
base_dir <- "~/Bureau/DataPhenolep/OutputRelernnCluster"

files <- list.files(
  path = base_dir,
  pattern = ".*\\.PREDICT\\.BSCORRECTED\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

for (f in files) {
  
  file_name <- basename(f)
  file_name <- gsub("\\.PREDICT\\.BSCORRECTED\\.txt$", "", file_name)
  
  out_data <- file.path(dirname(f), paste0(file_name, ".quantiles.tsv"))
  out_quant <- file.path(dirname(f), paste0(file_name, ".quantile_thresholds.tsv"))
  
  # ── Skip si déjà fait ──
  if (file.exists(out_data) && file.exists(out_quant)) {
    cat("Skipping (already done):", f, "\n")
    next
  }
  
  cat("Processing:", f, "\n")
  
  data <- read.table(f, header = TRUE)
  
  # Vérification
  if (!"recombRate" %in% colnames(data)) {
    cat("  -> recombRate column not found, skipped\n")
    next
  }
  
  # ── Conversion en cM/Mb ──
  data$recomb_cM_Mb <- data$recombRate * 1e8
  
  # ── Quantiles (seuils) ──
  quantile_values <- quantile(
    data$recomb_cM_Mb,
    probs = seq(0, 1, 0.01),
    na.rm = TRUE
  )
  
  quantile_values <- unique(quantile_values) ## ATTENTION, si plusieurs valeurs identiques entre les quantiles alors n < 100
  
  # ── Attribution quantile 1–100 ──
  data$recomb_quantile_100 <- cut(
    data$recomb_cM_Mb,
    breaks = quantile_values,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  # ── Export data ──
  write.table(
    data,
    file = out_data,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ── Export seuils ──
  write.table(
    quantile_values,
    file = out_quant,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
}