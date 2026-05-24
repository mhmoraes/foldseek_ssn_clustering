#!/usr/bin/env Rscript
library(tidyverse)
library(igraph)
library(tidygraph)

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- "results/foldseek_alignments.tsv"
out_dir   <- "results"

if (length(args) > 0) {
  for (i in 1:(length(args)-1)) {
    if (args[i] == "--input")  input_tsv <- args[i+1]
    if (args[i] == "--outdir") out_dir   <- args[i+1]
  }
}

# Load complete data matrix
alignment_columns <- c("query", "target", "lddt", "prob", "evalue", "fident", "alntmscore")
raw_data <- read_delim(input_tsv, delim = "\t", col_names = alignment_columns, show_col_types = FALSE)
all_nodes <- unique(c(raw_data$query, raw_data$target))

# --- STEP 1: CALCULATE ORIGINAL RUN METRICS ---
threshold_range <- seq(0.40, 0.95, by = 0.02)
original_modularity <- numeric(length(threshold_range))
unbiased_modularity  <- numeric(length(threshold_range))

idx <- 1
for (cutoff in threshold_range) {
  # Raw slice
  edges_raw <- raw_data %>% filter(lddt >= cutoff)
  if (nrow(edges_raw) > 0) {
    g_raw <- graph_from_data_frame(d = edges_raw %>% select(query, target, lddt), directed = FALSE, vertices = tibble(name = all_nodes))
    original_modularity[idx] <- modularity(cluster_louvain(g_raw))
  }
  
  # Unbiased slice: Drop edges that are purely driven by high sequence similarity (>80% identity)
  # This mathematically simulates how the network behaves without overrepresented homolog groups
  edges_unbiased <- raw_data %>% filter(lddt >= cutoff) %>% filter(fident < 0.80)
  if (nrow(edges_unbiased) > 0) {
    g_unbiased <- graph_from_data_frame(d = edges_unbiased %>% select(query, target, lddt), directed = FALSE, vertices = tibble(name = all_nodes))
    unbiased_modularity[idx] <- modularity(cluster_louvain(g_unbiased))
  }
  idx <- idx + 1
}

# --- STEP 2: EXTRACT OPTMIZED METRICS COMPARISONS ---
opt_idx_raw      <- admitting_max_idx <- which.max(original_modularity)
opt_cutoff_raw   <- threshold_range[opt_idx_raw]
max_mod_raw      <- original_modularity[opt_idx_raw]

opt_idx_unbiased <- which.max(unbiased_modularity)
opt_cutoff_unb   <- threshold_range[opt_idx_unbiased]
max_mod_unb      <- unbiased_modularity[opt_idx_unbiased]

# --- STEP 3: EXPORT TECHNICAL BIAS REPORT ---
report_path <- file.path(out_dir, "bias_analysis_report.txt")
report_conn <- file(report_path, "w")
writeLines("========================================================================", report_conn)
writeLines("STRUCTURAL SIMILARITY NETWORK SAMPLING BIAS EVALUATION REPORT", report_conn)
writeLines("========================================================================", report_conn)
writeLines(sprintf("Date evaluated: %s", Sys.time()), report_conn)
writeLines("", report_conn)
writeLines("1. RAW NETWORK PROPERTIES (INCLUDES ALL SEQUENCES):", report_conn)
writeLines(sprintf("   -> Mathematically Optimal lDDT Threshold: %s", opt_cutoff_raw), report_conn)
writeLines(sprintf("   -> Peak Network Modularity (M):          %s", round(max_mod_raw, 4)), report_conn)
writeLines("", report_conn)
writeLines("2. UNBIASED NETWORK PROPERTIES (THINNING OUT HIGH REDUNDANCY CORE VIA IDENTITY < 80%):", report_conn)
writeLines(sprintf("   -> Mathematically Optimal lDDT Threshold: %s", opt_cutoff_unb), report_conn)
writeLines(sprintf("   -> Peak Network Modularity (M):          %s", round(max_mod_unb, 4)), report_conn)
writeLines("", report_conn)
writeLines("3. CONTEXTUAL DIAGNOSTIC INTERPRETATION:", report_conn)

if (abs(opt_cutoff_raw - opt_cutoff_unb) <= 0.04) {
  writeLines("   [CONCLUSION] The optimal threshold is STABLE. Your dataset's structural clustering", report_conn)
  writeLines("   is robust and driven by genuine ancestral fold topology rather than database sampling bias.", report_conn)
} else {
  writeLines("   [WARNING] The optimal threshold shifted significantly. Your clustering runs are", report_conn)
  writeLines("   vulnerable to heavy organismal sampling artifacts. Consider utilizing the thinned dataset.", report_conn)
}
close(report_conn)

# --- STEP 4: GENERATE DIAGNOSTIC COMPARISON PLOT ---
plot_df <- tibble(
  Threshold = rep(threshold_range, 2),
  Modularity = c(original_modularity, unbiased_modularity),
  Dataset = c(rep("Raw Network (All Data)", length(threshold_range)), 
              rep("Unbiased Network (Sequence Identity < 80%)", length(threshold_range)))
)

comparison_chart <- ggplot(plot_df, aes(x = Threshold, y = Modularity, color = Dataset, linetype = Dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = opt_cutoff_raw, color = "darkblue", linetype = "dotted", alpha = 0.7) +
  geom_vline(xintercept = opt_cutoff_unb, color = "darkred", linetype = "dotted", alpha = 0.7) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Modularity Sensitivity Analysis: Sampling Bias Audit",
    x = "Foldseek lDDT Cutoff Threshold",
    y = "Network Modularity (M)",
    caption = "Dotted vertical lines indicate locked peak optimization coordinates."
  )

ggsave(file.path(out_dir, "bias_comparison_matrix.png"), plot = comparison_chart, width = 8, height = 6, dpi = 300)
cat(sprintf("[%s] [SUCCESS] Sensitivity audit reports generated cleanly.\n", Sys.time()))