# High-Throughput Automated Fault-Tolerant Foldseek SSN Pipeline

A pipeline managed via Snakemake that processes folders of protein PDB coordinate files, executes local all-vs-all structural alignments, optimizes subgroup partitioning by maximizing network modularity ($M$), and outputs publication-grade network maps and sequence identity histograms.

## 1. Installation & Execution Strategy

Ensure that Mamba or Conda is active on your Linux/macOS workstation.

```bash
# Install Snakemake inside your base profile if not already available
conda install -c conda-forge snakemake

# Run the complete pipeline automatically handling environment builds on the fly
snakemake --use-conda --cores 4
