Here is the completely updated, production-grade `README.md` file for your repository. It integrates the newly implemented **Down-sampling Sensitivity Analysis** and documents every single deliverable inside your `results/` folder so your workspace remains fully self-explanatory.

Overwrite your current `README.md` file with this comprehensive layout:

#### `README.md`

```markdown
# High-Throughput Automated Foldseek SSN & Sampling Bias Pipeline

An automated, reproducible, data-driven bioinformatics platform managed via Snakemake. This pipeline ingests unaligned protein PDB coordinate assets, executes high-speed all-vs-all structural alignments, optimizes subgroup partitioning by maximizing global network modularity ($M$), and evaluates whether database overrepresentation (sampling bias) is distorting the resulting network layout.

---

## 1. Pipeline Architecture & Execution Flow

The workflow is split into three decoupled, environment-isolated processing rules:

1. **`run_foldseek_alignment` (Python Backend):** Indexes raw PDB structures, runs high-sensitivity structural alignments, and generates a square lDDT similarity matrix.
2. **`optimize_and_render_ssn` (R Graphics Engine):** Sweeps through lDDT thresholds to find the peak modularity coordinate, partitions proteins into optimized subfamilies, and renders publication-ready vector graph maps.
3. **`evaluate_clustering_bias` (R Sensitivity Audit):** Simulates a sequence redundancy filter (identity < 80%) to check if overrepresented organisms are skewing the structural network architecture.

---

## 2. Installation & Quick Start

Ensure that Mamba or Conda is installed and active on your macOS/Linux workstation.

```bash
# Clone the repository and navigate to your project root folder
cd structure_clustering

# Configure strict channel priorities to guarantee rapid environment builds
conda config --set channel_priority strict

# Launch the complete automated pipeline end-to-end
snakemake --use-conda --cores 4

```

*Note: If your local runtime engine encounters a pipeline folder lock, clear it by executing `snakemake --unlock` before restarting.*

---

## 3. Input Data Specifications

* **Directory Target:** `extracted_domains/`
* **File Requirements:** A folder populated by individual, unaligned structural coordinate files in standard `.pdb` or `.cif` formats.
* **Biological Context:** Optimized for evaluating evolutionary structural drift across large enzyme superfamilies (e.g., Deaminases) deep within the sequence "Twilight Zone."

---

## 4. Exhaustive Deliverables Ledger (`results/`)

Every file written into the centralized `results/` folder is designed to be fully complementary:

| File Name | Format | Primary Computational & Scientific Utility |
| --- | --- | --- |
| `foldseek_alignments.tsv` | Tabular Edge List | The raw, unrestricted pairwise alignment output fields (`query`, `target`, `lddt`, `prob`, `evalue`, `fident`, `alntmscore`). Serves as the master database for downstream filters. |
| `structural_similarity_matrix.tsv` | Symmetric Tensor | A pristine, square distance matrix populated by mutual lDDT match scores. Ready for direct ingestion into machine learning frameworks (PCA, UMAP, t-SNE, or Scikit-Learn classifiers). |
| `optimized_subfamily_assignments.tsv` | Flat Spreadsheet | The master mapping ledger linking each protein node to its optimized cluster ID and its internal connection hubness (`node_degree`). Sorting by degree isolates core structural templates for lab testing. |
| `optimized_subfamilies_network_map.pdf` | Vector Graphic | A high-resolution network map colored by subfamily and sized by hubness. Saved as geometric paths ready for direct modifications inside Adobe Illustrator or Inkscape. |
| `sequence_identity_distribution.png` | Raster Chart | A multimodal frequency histogram plotting pairwise fractional identity (`fident`). Distinct peaks below 20% confirm structural conservation despite extensive primary sequence drift. |
| `ssn_graph_object.rds` | Serialized Binary | A live, compressed R data list containing your network model, layout coordinates, and optimal cutoff. Allows you to load your workspace back into RStudio instantly to tweak graphics without re-running scripts. |
| `bias_comparison_matrix.png` | Diagnostic Line | **[NEW]** A dual-line plot tracing raw network modularity against an unbiased network (identity < 80%). Overlapping trajectories prove your structural layout is robust against database sampling errors. |
| `bias_analysis_report.txt` | Text Document | **[NEW]** A formal mathematical report evaluating the stability of your threshold under sequence thinning, providing an automated conclusion on whether bias is skewing your results. |

---

## 5. Interactive Sandbox: Tweaking Graphics Natively

To modify your visualization styling, layouts, or color schemes dynamically without re-running the heavy computational backend steps, open an interactive RStudio session from your root directory and leverage the serialized `.rds` data package:

```R
library(tidyverse)
library(ggraph)
library(colorspace)

# Instantly restore the live pipeline network data structure
ssn_data <- readRDS("results/ssn_graph_object.rds")

# Re-plot using an alternative palette (e.g., Viridis)
ggraph(ssn_data$graph, layout = "manual", x = ssn_data$coordinates[,1], y = ssn_data$coordinates[,2]) +
  geom_edge_diagonal0(color = "grey30", alpha = 0.4, show.legend = FALSE) +
  geom_node_point(aes(color = community_id, size = node_degree), alpha = 0.9) +
  scale_color_manual(values = colorspace::sequential_hcl(length(unique(V(ssn_data$graph)$community_id)), "Viridis")) +
  theme_graph()

```

```

---

### Push the Updated Repository to GitHub

Since your local repository is fully up to date, stage your modified files and stream them to your GitHub cloud repository:

```bash
# 1. Verify that your modified files are tracked cleanly
git status

# 2. Stage your updated README and your new sensitivity script
git add README.md scripts/evaluate_clustering_bias.R Snakefile

# 3. Commit using clean engineering tracking tokens
git commit -m "docs: upgrade README to explicitly detail the down-sampling clustering bias sensitivity deliverables"

# 4. Push up to your remote master timeline
git push origin main

```