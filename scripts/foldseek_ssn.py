import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import networkx as nx

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

def run_local_foldseek(input_dir: Path, output_dir: Path, cores: int) -> Path:
    checkpoint_db_dir = output_dir / "indexed_database"
    checkpoint_db_dir.mkdir(parents=True, exist_ok=True)
    
    local_db = checkpoint_db_dir / "local_protein_db"
    raw_alignments = checkpoint_db_dir / "alignment_results"
    output_tsv = output_dir / "foldseek_alignments.tsv"

    if not (local_db.exists() and Path(str(local_db) + ".index").exists()):
        subprocess.run(["foldseek", "createdb", str(input_dir), str(local_db)], check=True, stdout=subprocess.PIPE)
    
    if not (raw_alignments.exists() and Path(str(raw_alignments) + ".index").exists()):
        subprocess.run(["foldseek", "search", str(local_db), str(local_db), str(raw_alignments), str(checkpoint_db_dir), "-s", "9.5", "--threads", str(cores), "--alignment-type", "1", "--prefilter-mode", "1", "-a"], check=True, stdout=subprocess.PIPE)
    
    if not (output_tsv.exists() and output_tsv.stat().st_size > 0):
        format_fields = "query,target,lddt,prob,evalue,fident,alntmscore"
        subprocess.run(["foldseek", "convertalis", str(local_db), str(local_db), str(raw_alignments), str(output_tsv), "--format-output", format_fields], check=True, stdout=subprocess.PIPE)
        
    return output_tsv

def generate_matrix(foldseek_tsv: Path, matrix_tsv: Path) -> None:
    column_names = ["query", "target", "lddt", "prob", "evalue", "fident", "alntmscore"]
    df = pd.read_csv(foldseek_tsv, sep="\t", header=None, names=column_names)
    nodes = sorted(list(set(df["query"].unique()).union(set(df["target"].unique()))))
    
    matrix = pd.DataFrame(np.eye(len(nodes)), index=nodes, columns=nodes)
    for _, row in df.iterrows():
        matrix.loc[row["query"], row["target"]] = float(row["lddt"])
        matrix.loc[row["target"], row["query"]] = float(row["lddt"])
        
    matrix.to_csv(matrix_tsv, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-c", "--cores", type=int, default=4)
    args = parser.parse_args()
    
    out_path = Path(args.output_dir)
    align_tsv = run_local_foldseek(Path(args.input_dir), out_path, args.cores)
    generate_matrix(align_tsv, out_path / "structural_similarity_matrix.tsv")
