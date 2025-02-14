#!/usr/bin/env python3
"""
This script iterates through all files in the specified folder that match the pattern:
    Phe_*_pval_up_to_0.1.gz
For each file, it:
  - Reads variant data from the Tabix-indexed file (and its corresponding .tbi file).
  - Computes minor allele frequency (MAF) and -log10(p) values.
  - Generates a QQ plot (observed vs. expected -log10(p)) grouped by MAF bins.
  - Annotates the plot with Genomic Control lambda values.
  - Saves the QQ plot as a PNG image with the same base name as the input file.
  
Each input file is assumed to have a corresponding Tabix index file named by appending ".tbi".
"""

import os
import glob
import pysam
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.stats import beta
from scipy.stats import rankdata

# ----------------------------------------------------------------------------
# 1) CONFIGURATION
# ----------------------------------------------------------------------------

# Directory containing your .gz files
base_dir = "/Users/hritvik/genomics-backend/DATABASE"  # Adjust as needed

# We use the same directory for output images
output_dir = base_dir

# File pattern: only files that have threshold 0.1 in their name.
file_pattern = os.path.join(base_dir, "Phe_*_pval_up_to_0.1.gz")

# Process only autosomal chromosomes 1-22
CHROMOSOMES = [str(i) for i in range(1, 23)]

# Define MAF categories along with colors for the QQ plot
MAF_RANGE_COLORS = {
    "0 ≤ MAF < 1e-4": {
        "range": (0, 0.0001),
        "color": "#ff7f0e"
    },
    "1e-4 ≤ MAF < 5e-4": {
        "range": (0.0001, 0.0005),
        "color": "#ffa500"
    },
    "5e-4 ≤ MAF < 0.02": {
        "range": (0.0005, 0.02),
        "color": "#a9a9a9"
    },
    "0.02 ≤ MAF < 0.50": {
        "range": (0.02, 0.5),
        "color": "#4169e1"
    }
}

# ----------------------------------------------------------------------------
# 2) HELPER FUNCTIONS
# ----------------------------------------------------------------------------

def fetch_tabix_data(chrom, gz_path, tbi_path):
    """
    Fetch variant records for a given chromosome from the Tabix file.
    Parses the following columns (0-indexed):
      - position: column 2
      - p-value: column 14
      - allele frequency (af): column 6
    Returns a list of dictionaries with keys: chrom, pos, pval, af.
    """
    results = []
    try:
        with pysam.TabixFile(gz_path, index=tbi_path) as tbx:
            for row in tbx.fetch(chrom):
                fields = row.strip().split('\t')
                if len(fields) < 19:
                    continue
                try:
                    pos_str  = fields[2]
                    pval_str = fields[14]
                    af_str   = fields[6]
                    
                    # Skip rows with missing values
                    if pval_str in ("NA", "") or af_str in ("NA", ""):
                        continue
                    
                    pos  = int(float(pos_str))
                    pval = float(pval_str)
                    af   = float(af_str)
                    
                    # Only accept variants with valid p-value and allele frequency
                    if pval > 0 and 0 < af < 1:
                        results.append({
                            "chrom": chrom,
                            "pos": pos,
                            "pval": pval,
                            "af": af
                        })
                except ValueError:
                    continue
    except Exception as e:
        print(f"Warning: Could not fetch data for chromosome {chrom} from {gz_path}.\n{e}")
    return results

def fetch_all_variants(gz_path, tbi_path):
    """
    Fetch variant data from all chromosomes and return as a pandas DataFrame.
    """
    all_data = []
    for c in CHROMOSOMES:
        chunk = fetch_tabix_data(c, gz_path, tbi_path)
        all_data.extend(chunk)
    return pd.DataFrame(all_data)

def calculate_maf(af):
    """
    Given allele frequency (af), return the minor allele frequency (MAF).
    """
    return min(af, 1 - af)

def calculate_gc_lambda(pvals, percentile=0.5):
    """
    Calculate the Genomic Control lambda at a specified percentile.
    pvals: 1D array-like of p-values.
    """
    if len(pvals) == 0:
        return float("nan")
    sorted_pvals = np.sort(pvals)
    obs = -np.log10(sorted_pvals)
    exp = -np.log10(np.arange(1, len(pvals) + 1) / (len(pvals) + 1))
    idx = int(len(pvals) * percentile)
    if idx < 1 or idx >= len(pvals):
        return float("nan")
    return obs[idx] / exp[idx]

def qq_unif_trace(pvals, label=None, color='blue',
                  already_transformed=False, should_thin=True,
                  thin_obs_places=2, thin_exp_places=2,
                  draw_conf=False, conf_points=1000, conf_alpha=0.05,
                  marker_size=3, marker_opacity=0.7):
    """
    Generate Plotly trace(s) for a QQ plot based on p-values.
    
    The expected values are computed as:
         expected = -log10( (rank - 0.5) / (N + 1) )
    where N is the total number of p-values.
    
    Optionally returns a confidence envelope trace.
    """
    pvals = np.array(pvals)
    N = len(pvals)
    if N == 0:
        raise ValueError("p-value vector is empty")
    
    if not already_transformed:
        if np.any(pvals == 0):
            raise ValueError("p-value vector contains zeros")
        obs = -np.log10(pvals)
        ranks = rankdata(pvals, method='ordinal')
        exp = -np.log10((ranks - 0.5) / (N + 1))
    else:
        obs = pvals
        ranks = rankdata(pvals, method='ordinal')
        exp = -np.log10(((N + 1) - ranks - 0.5) / (N + 1))
    
    if should_thin:
        df_temp = pd.DataFrame({
            'exp': np.round(exp, thin_exp_places),
            'obs': np.round(obs, thin_obs_places)
        })
        df_temp = df_temp.drop_duplicates()
        exp = df_temp['exp'].values
        obs = df_temp['obs'].values
    
    scatter_trace = go.Scatter(
        x=exp,
        y=obs,
        mode='markers',
        marker=dict(color=color, size=marker_size, opacity=marker_opacity),
        name=label if label is not None else ""
    )
    
    traces = []
    if draw_conf:
        points = min(conf_points, N)
        x_upper = []
        y_upper = []
        for i in range(1, points + 1):
            x_val = -np.log10((i - 0.5) / (N + 1))
            upper = -np.log10(beta.ppf(1 - conf_alpha / 2, i, N - i + 1))
            x_upper.append(x_val)
            y_upper.append(upper)
        x_lower = x_upper[::-1]
        y_lower = []
        for i in range(1, points + 1):
            lower = -np.log10(beta.ppf(conf_alpha / 2, i, N - i + 1))
            y_lower.append(lower)
        y_lower = y_lower[::-1]
        x_ci = x_upper + x_lower
        y_ci = y_upper + y_lower
        ci_trace = go.Scatter(
            x=x_ci,
            y=y_ci,
            fill='toself',
            fillcolor='lightgray',
            line=dict(color='lightgray'),
            hoverinfo='skip',
            showlegend=False
        )
        traces.append(ci_trace)
    
    traces.append(scatter_trace)
    return traces

# ----------------------------------------------------------------------------
# 3) MAIN PROCESSING LOOP
# ----------------------------------------------------------------------------

def main():
    # Find all files matching the pattern (with pval_up_to_0.1)
    file_list = glob.glob(file_pattern)
    if not file_list:
        print("No files matching the pattern were found in the folder.")
        return

    print(f"Found {len(file_list)} file(s). Processing each one...")
    
    for gz_file_path in file_list:
        tbi_file_path = gz_file_path + ".tbi"  # Assuming the index file is gz_file.tbi
        base_name = os.path.basename(gz_file_path)  # e.g., Phe_401_1.AFR.gwama_pval_up_to_0.1.gz
        out_base = base_name.replace(".gz", "")      # e.g., Phe_401_1.AFR.gwama_pval_up_to_0.1
        output_file = os.path.join(output_dir, out_base + ".png")
        
        print(f"\nProcessing file: {base_name}")
        if not os.path.exists(tbi_file_path):
            print(f"  Warning: Tabix index file {tbi_file_path} not found. Skipping.")
            continue
        
        print("  Reading variant data ...")
        df = fetch_all_variants(gz_file_path, tbi_file_path)
        print(f"  Fetched {len(df):,} variants.")
        
        if df.empty:
            print("  No valid data found. Skipping plot generation for this file.")
            continue
        
        # Compute MAF and -log10(p)
        df["maf"] = df["af"].apply(calculate_maf)
        df["logp"] = -np.log10(df["pval"])
        
        # Calculate overall GC lambda values for annotation
        pvals_all = df["pval"].values
        lambdas = {
            "0.5": calculate_gc_lambda(pvals_all, 0.5),
            "0.1": calculate_gc_lambda(pvals_all, 0.1),
            "0.01": calculate_gc_lambda(pvals_all, 0.01),
            "0.001": calculate_gc_lambda(pvals_all, 0.001),
        }
        
        # Prepare QQ plot traces
        plot_traces = []
        sorted_all = np.sort(pvals_all)
        N_all = len(sorted_all)
        obs_all = -np.log10(sorted_all)
        exp_all = -np.log10((np.arange(1, N_all + 1) - 0.5) / (N_all + 1))
        max_val = max(np.max(obs_all), np.max(exp_all))
        
        # Add the identity line (y = x)
        identity_trace = go.Scatter(
            x=[0, max_val],
            y=[0, max_val],
            mode="lines",
            line=dict(color="black", dash="dash", width=1),
            showlegend=False,
            hoverinfo="skip"
        )
        plot_traces.append(identity_trace)
        
        # Create QQ traces for each MAF category
        for label, info in MAF_RANGE_COLORS.items():
            min_maf, max_maf = info["range"]
            color = info["color"]
            cat_df = df[(df["maf"] >= min_maf) & (df["maf"] < max_maf)]
            if cat_df.empty:
                continue
            
            sorted_cat = np.sort(cat_df["pval"].values)
            traces = qq_unif_trace(sorted_cat,
                                   label=f"{label} ({len(cat_df):,})",
                                   color=color,
                                   draw_conf=False,
                                   marker_size=3,
                                   marker_opacity=0.7)
            plot_traces.extend(traces)
        
        # Build GC lambda annotation text
        lambda_text = (
            f"GC lambda 0.5: {lambdas['0.5']:.3f}<br>"
            f"GC lambda 0.1: {lambdas['0.1']:.3f}<br>"
            f"GC lambda 0.01: {lambdas['0.01']:.3f}<br>"
            f"GC lambda 0.001: {lambdas['0.001']:.3f}<br>"
            "(Genomic Control lambda at specified percentiles)"
        )
        
        # Create the Plotly layout
        layout = go.Layout(
            title=f"QQ Plot for {base_name}",
            xaxis=dict(
                title="Expected -log10(p)",
                range=[0, max_val],
                showline=True,
                linecolor="black",
                mirror=True,
                gridcolor="#E2E2E2"
            ),
            yaxis=dict(
                title="Observed -log10(p)",
                range=[0, max_val],
                showline=True,
                linecolor="black",
                mirror=True,
                gridcolor="#E2E2E2"
            ),
            legend=dict(
                x=0.03,
                y=0.97,
                bgcolor="rgba(255,255,255,0.8)"
            ),
            annotations=[
                dict(
                    x=0.03,
                    y=-0.15,
                    xref="paper",
                    yref="paper",
                    text=lambda_text,
                    showarrow=False
                )
            ],
            width=700,
            height=700,
            plot_bgcolor="white"
        )
        
        fig = go.Figure(data=plot_traces, layout=layout)
        
        try:
            fig.write_image(output_file, scale=2)
            print(f"  QQ plot saved to: {output_file}")
        except Exception as e:
            print(f"  Error saving QQ plot for {base_name}:\n{e}")

if __name__ == "__main__":
    main()
