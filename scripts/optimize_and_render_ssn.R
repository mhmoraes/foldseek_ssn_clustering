#!/usr/bin/env Rscript
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(colorspace)

# Parse parameters from the Snakemake context
args <- commandArgs(trailingOnly = TRUE)
input_tsv <- "results/foldseek_alignments.tsv"
out_dir   <- "results"

if (length(args) > 0) {
  for (i in 1:(length(args)-1)) {
    if (args[i] == "--input")  input_tsv <- args[i+1]
    if (args[i] == "--outdir") out_dir   <- args[i+1]
  }
}

# Ingest data schema
alignment_columns <- c("query", "target", "lddt", "prob", "evalue", "fident", "alntmscore")
raw_data <- read_delim(input_tsv, delim = "\t", col_names = alignment_columns, show_col_types = FALSE)
all_nodes <- unique(c(raw_data$query, raw_data$target))

# --- STEP 2A: SWEEP PARAMETER CONTUNUUM ---
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

# ==============================================================================
# --- STEP 2B: NETWORK TRANSFORMATION (SELF-LOOPS PURGED) ---
# ==============================================================================
# 1. Filter edges strictly by your mathematically optimized cutoff,
# AND explicitly drop all self-loops where a protein aligns to itself!
optimized_edges <- raw_data %>% 
  filter(lddt >= optimal_cutoff) %>%
  filter(query != target)  # CRITICAL FIX: Eliminates ghost degrees from self-matching

# 2. Build the final igraph structure ensuring ALL 870 original nodes exist
final_igraph <- graph_from_data_frame(
  d = optimized_edges %>% select(query, target, lddt), 
  directed = FALSE,
  vertices = tibble(name = all_nodes)
)

# 3. Calculate degrees using igraph's native function on the clean graph
true_degrees <- degree(final_igraph)

# 4. Convert safely to tidygraph with absolute in-place attribute mapping
final_tidygraph <- as_tbl_graph(final_igraph) %>%
  activate(nodes) %>%
  mutate(community_id = as_factor(group_louvain())) %>%
  mutate(node_degree = as.integer(true_degrees[name]))

# 5. Generate the layout coordinates safely using the stress engine
final_layout <- create_layout(final_tidygraph, layout = "stress", bbox = 15)

# 6. Extract and export the clean spreadsheet ledger directly from the layout
node_ledger <- as_tibble(final_layout) %>%
  rename(protein_id = name) %>%
  select(protein_id, node_degree, community_id)

write_tsv(node_ledger, file.path(out_dir, "optimized_subfamily_assignments.tsv"))

# ==============================================================================
# --- STEP 2C: SERIALIZED RDS DATA BACKUP OBJECT ---
# ==============================================================================
graph_package <- list(
  graph       = final_tidygraph,
  layout      = final_layout,
  cutoff      = optimal_cutoff
)
saveRDS(graph_package, file = file.path(out_dir, "ssn_graph_object.rds"))
cat(sprintf("[%s] [SUCCESS] Exported synchronized R data package to results/ssn_graph_object.rds\n", Sys.time()))

# ==============================================================================
# --- STEP 2D: RENDER & SAVE PUBLICATION GRAPHICS ---
# ==============================================================================
num_clades <- length(unique(node_ledger$community_id))

ssn_plot <- ggraph(final_layout) +
  geom_edge_diagonal0(aes(alpha = lddt), color = "grey40", show.legend = FALSE) +
  geom_node_point(aes(color = community_id, size = node_degree), alpha = 0.85) +
  scale_color_manual(values = colorspace::rainbow_hcl(num_clades, c = 75, l = 60)) +
  scale_size_continuous(name = "Hubness (Node Degree)", range = c(1.5, 7.5)) +
  scale_edge_alpha_continuous(range = c(0.2, 0.9)) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold")) +
  guides(color = "none") +
  labs(
    title = "Optimized Structural Similarity Network (SSN)", 
    subtitle = sprintf("Foldseek Engine | Cutoff Maxima lDDT >= %s | Subfamilies Resolved: %d", optimal_cutoff, num_clades)
  )

ggsave(file.path(out_dir, "optimized_subfamilies_network_map.pdf"), plot = ssn_plot, width = 11, height = 8.5, device = "pdf")

# ==============================================================================
# --- STEP 2E: SEQUENCE IDENTITY HISTOGRAM ---
# ==============================================================================
png(file.path(out_dir, "sequence_identity_distribution.png"), width = 800, height = 600, res = 120)
hist(raw_data$fident, breaks = 50, main = "Sequence Identity Distribution", xlab = "Fractional Identity (fident)", col = "skyblue", border = "black")
dev.off()