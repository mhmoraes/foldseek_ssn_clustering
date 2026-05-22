#!/usr/bin/env Rscript
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(colorspace)

# Parse standard bash-orchestration string flags from the Snakemake context
args <- commandArgs(trailingOnly = TRUE)
input_tsv  <- args[which(args == "--input") + 1]
out_dir    <- args[which(args == "--outdir") + 1]

# Ingest data schema
alignment_columns <- c("query", "target", "lddt", "prob", "evalue", "fident", "alntmscore")
raw_data <- read_delim(input_tsv, delim = "\t", col_names = alignment_columns, show_col_types = FALSE)
all_nodes <- unique(c(raw_data$query, raw_data$target))

# --- STEP 2A: AUTOMATED PARALLEL PARAMETER SENSITIVITY CONTUNUUM SWEEP ---
threshold_range <- seq(0.40, 0.95, by = 0.02)
best_modularity <- -1
optimal_cutoff  <- 0.70

for (cutoff in threshold_range) {
  slice_edges <- raw_data %>% filter(lddt >= cutoff)
  if (nrow(slice_edges) == 0) next
  slice_graph <- graph_from_data_frame(d = slice_edges %>% select(query, target, lddt), directed = FALSE, vertices = tibble(name = all_nodes))
  mod_score <- modularity(cluster_louvain(slice_graph))
  if (mod_score > best_modularity) {
    best_modularity <- mod_score
    optimal_cutoff  <- cutoff
  }
}

# --- STEP 2B: NETWORK PARTITION MANIFEST EXTRACTION ---
optimized_edges <- raw_data %>% filter(lddt >= optimal_cutoff)
final_igraph <- graph_from_data_frame(d = optimized_edges %>% select(query, target, lddt), directed = FALSE, vertices = tibble(name = all_nodes))
coords <- layout_components(final_igraph, layout = layout_with_stress)

final_tidygraph <- as_tbl_graph(final_igraph) %>%
  activate(nodes) %>%
  mutate(node_degree = centrality_degree(), community_id = as_factor(group_louvain()))

# Export clean spreadsheet assignments list
node_ledger <- final_tidygraph %>% activate(nodes) %>% as_tibble() %>% rename(protein_id = name)
write_tsv(node_ledger, file.path(out_dir, "optimized_subfamily_assignments.tsv"))

# --- STEP 2C: HIGH-CONTRAST COMPONENT CANVAS RENDERING ---
num_clades <- length(unique(node_ledger$community_id))
ssn_plot <- ggraph(final_tidygraph, layout = "manual", x = coords[,1], y = coords[,2]) +
  geom_edge_diagonal0(aes(alpha = lddt), color = "grey40", show.legend = FALSE) +
  geom_node_point(aes(color = community_id, size = node_degree), alpha = 0.85) +
  scale_color_manual(values = colorspace::rainbow_hcl(num_clades, c = 75, l = 60)) +
  scale_size_continuous(name = "Hubness (Node Degree)", range = c(1.5, 7.5)) +
  scale_edge_alpha_continuous(range = c(0.2, 0.9)) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold")) +
  guides(color = "none") +
  labs(title = "Optimized Structural Similarity Network (SSN)", subtitle = sprintf("Foldseek Adjacency Engine | Cutoff Maxima lDDT >= %s | Subfamilies Resolved: %d", optimal_cutoff, num_clades))

ggsave(file.path(out_dir, "optimized_subfamilies_network_map.pdf"), plot = ssn_plot, width = 11, height = 8.5, device = "pdf")

# --- STEP 2D: MULTI-MODAL DISTRIBUTIONS SUMMARY HISTOGRAM ---
png(file.path(out_dir, "sequence_identity_distribution.png"), width = 800, height = 600, res = 120)
hist(raw_data$fident, breaks = 50, main = "Sequence Identity Distribution", xlab = "Fractional Identity (fident)", col = "skyblue", border = "black")
dev.off()
