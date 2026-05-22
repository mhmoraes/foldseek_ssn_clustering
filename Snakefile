"""Snakemake pipeline orchestrating data-driven Foldseek SSN generation and R modularity graphics."""

import os

# Define environment variables or configuration defaults
INPUT_PDB_DIR = "extracted_domains"
OUTPUT_DIR = "results"
ALLOCATED_CORES = 4

rule all:
    """Target rule specifying the final structural analytics deliverables."""
    input:
        matrix = os.path.join(OUTPUT_DIR, "structural_similarity_matrix.tsv"),
        assignments = os.path.join(OUTPUT_DIR, "optimized_subfamily_assignments.tsv"),
        network_map = os.path.join(OUTPUT_DIR, "optimized_subfamilies_network_map.pdf"),
        histogram = os.path.join(OUTPUT_DIR, "sequence_identity_distribution.png"),
        rds_object = os.path.join(OUTPUT_DIR, "ssn_graph_object.rds") # ADDED: Live R Graphic Sandbox File


rule run_foldseek_alignment:
    """Executes high-speed all-vs-all structural alignments utilizing local Foldseek databases."""
    input:
        pdb_dir = INPUT_PDB_DIR
    output:
        alignments = os.path.join(OUTPUT_DIR, "foldseek_alignments.tsv"),
        matrix = os.path.join(OUTPUT_DIR, "structural_similarity_matrix.tsv")
    params:
        script = "scripts/foldseek_ssn.py",
        cores = ALLOCATED_CORES
    conda:
        "envs/foldseek_env.yml"
    log:
        os.path.join(OUTPUT_DIR, "logs/foldseek_alignment.log")
    shell:
        """
        python {params.script} \
          -i {input.pdb_dir} \
          -o {OUTPUT_DIR} \
          -c {params.cores} > {log} 2>&1
        """


rule optimize_and_render_ssn:
    """Sweeps modularity values, extracts optimized subfamilies, and generates publication plots."""
    input:
        alignments = os.path.join(OUTPUT_DIR, "foldseek_alignments.tsv")
    output:
        assignments = os.path.join(OUTPUT_DIR, "optimized_subfamily_assignments.tsv"),
        network_map = os.path.join(OUTPUT_DIR, "optimized_subfamilies_network_map.pdf"),
        histogram = os.path.join(OUTPUT_DIR, "sequence_identity_distribution.png"),
        rds_object = os.path.join(OUTPUT_DIR, "ssn_graph_object.rds") # ADDED: Passthrough destination flag
    params:
        script = "scripts/optimize_and_render_ssn.R"
    conda:
        "envs/r_graphics_env.yml"
    log:
        os.path.join(OUTPUT_DIR, "logs/r_graphics_analysis.log")
    shell:
        """
        Rscript {params.script} \
          --input {input.alignments} \
          --outdir {OUTPUT_DIR} > {log} 2>&1
        """