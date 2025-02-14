# import pandas as pd
# import numpy as np
# import plotly.graph_objects as go
# import pysam
# from pathlib import Path
# import glob

# # Constants
# CHROMOSOMES = [str(i) for i in range(1, 23)]
# ALTERNATING_COLORS = ['#DC2626', '#2563EB']  # Dark red and blue
# DOWNSAMPLE_LIMIT = 50000

# def fetch_tabix_data(chrom, gz_file_path, tbi_file_path):
#     """Fetch data for a single chromosome from a tabix file."""
#     results = []
#     with pysam.TabixFile(gz_file_path, index=tbi_file_path) as tabix_file:
#         for row in tabix_file.fetch(str(chrom)):
#             fields = row.strip().split('\t')
#             try:
#                 pos = int(float(fields[2]))     # POS column (unchanged)
#                 pval = float(fields[7]) if fields[7] != 'NA' else None  # Changed from 14 to 7 for P column
#                 if pval is not None:
#                     results.append({
#                         'chrom': str(chrom),
#                         'pos': pos,
#                         'pval': pval,
#                     })
#             except (ValueError, IndexError) as e:
#                 print(f"Error processing row for chromosome {chrom}: {str(e)}")
#                 continue
#     return results

# def fetch_all_chromosomes(gz_file_path, tbi_file_path):
#     """Fetch data for all chromosomes from a single file pair."""
#     all_data = []
#     for chrom in range(1, 23):
#         all_data.extend(fetch_tabix_data(chrom, gz_file_path, tbi_file_path))
#     return pd.DataFrame(all_data)

# def process_dataframe(df):
#     """Process the dataframe and calculate manhattan plot coordinates."""
#     df['log_p'] = -np.log10(df['pval'])
#     df['chrom'] = df['chrom'].astype(int)
#     df = df.sort_values(['chrom', 'pos'])

#     # Calculate chromosome offsets with larger spacing
#     chrom_sizes = df.groupby('chrom')['pos'].max()
#     spacing = chrom_sizes.max() * 0.08
#     chrom_offsets = pd.Series(0, index=range(1, 23))
#     current_offset = 0

#     for chrom in range(1, 23):
#         chrom_offsets[chrom] = current_offset
#         if chrom in chrom_sizes.index:
#             current_offset += chrom_sizes[chrom] + spacing

#     # Apply offsets
#     df['x_manhattan'] = df.apply(lambda row: row['pos'] + chrom_offsets[row['chrom']], axis=1)

#     return df, chrom_sizes, chrom_offsets

# def create_manhattan_plot(df, chrom_sizes, chrom_offsets, output_path, height=150):
#     """Create and save a manhattan plot."""
#     # Calculate tick positions
#     tickvals = []
#     ticktext = []
#     for chrom in range(1, 23):
#         if chrom in chrom_sizes.index:
#             center = chrom_offsets[chrom] + chrom_sizes[chrom] / 2
#             tickvals.append(center)
#             ticktext.append(str(chrom))

#     # Create plot data
#     plot_data = []
#     for chrom in range(1, 23):
#         chrom_data = df[df['chrom'] == chrom]
#         color_index = (chrom - 1) % 2
#         plot_data.append(go.Scattergl(
#             x=chrom_data['x_manhattan'],
#             y=chrom_data['log_p'],
#             mode='markers',
#             marker=dict(
#                 color=ALTERNATING_COLORS[color_index],
#                 size=5,
#                 opacity=0.8
#             ),
#             hoverinfo='text',
#             text=chrom_data.apply(
#                 lambda row: f"Chromosome: {row['chrom']}<br>Position: {row['pos']}<br>-log10 p-value: {row['log_p']:.2f}<br>P-value: {10**-row['log_p']:.2e}",
#                 axis=1
#             ),
#             showlegend=False
#         ))

#     # Create layout
#     layout = go.Layout(
#         showlegend=False,
#         xaxis=dict(
#             title=None,
#             tickmode='array',
#             tickvals=tickvals,
#             ticktext=ticktext,
#             showline=False,
#             zeroline=False,
#             showgrid=False,
#             showticklabels=False,
#             ticks='',
#             ticklen=0,
#             tickfont=dict(size=10),
#         ),
#         yaxis=dict(
#             title=None,
#             range=[0, 10],
#             showline=False,
#             zeroline=False,
#             showgrid=False,
#             showticklabels=False,
#             ticks='',
#             ticklen=0,
#         ),
#         width=1500,
#         height=height,
#         plot_bgcolor='white',
#         paper_bgcolor='white',
#         margin=dict(
#             l=0,
#             r=0,
#             t=10,
#             b=30
#         ),
#     )

#     # Create and save figure
#     fig = go.Figure(data=plot_data, layout=layout)
#     fig.write_image(
#         output_path,
#         scale=3,
#         width=1500,
#         height=height
#     )
#     return fig

# def main():
#     # Directory containing the .gz files
#     data_dir = "/nfs/platlas_stor/tabix/"
#     output_dir = "/nfs/platlas_stor/manhatton_plot"
    
#     # Find all .gz files
#     gz_files = glob.glob(f"{data_dir}/*.gz")
#     # Filter out .tbi files from the list
#     gz_files = [f for f in gz_files if not f.endswith('.tbi')]
    
#     # Process each file pair
#     for gz_file_path in gz_files:
#         tbi_file_path = f"{gz_file_path}.tbi"
        
#         # Skip if .tbi file doesn't exist
#         if not Path(tbi_file_path).exists():
#             print(f"Warning: No .tbi file found for {gz_file_path}")
#             continue
            
#         # Generate output filename
#         file_name = Path(gz_file_path).stem
#         output_path = Path(output_dir) / f"manhattan_{file_name}.png"
        
#         print(f"Processing {file_name}...")
        
#         # Fetch and process data
#         df = fetch_all_chromosomes(gz_file_path, tbi_file_path)
#         df, chrom_sizes, chrom_offsets = process_dataframe(df)
        
#         # Create plot
#         create_manhattan_plot(df, chrom_sizes, chrom_offsets, output_path)
#         print(f"Saved plot to {output_path}")

# if __name__ == "__main__":
#     main()

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import pysam
from pathlib import Path
import glob

# Constants
CHROMOSOMES = [str(i) for i in range(1, 23)]
ALTERNATING_COLORS = ['#DC2626', '#2563EB']  # Dark red and blue
SIGNIFICANCE_THRESHOLD = 5e-8  # Genome-wide significance threshold

def process_dataframe(df):
    """Process the dataframe and calculate manhattan plot coordinates."""
    # Calculate -log10(p-value)
    df['log_p'] = -np.log10(df['pval'].clip(1e-300))  # Clip to prevent infinity
    df['chrom'] = df['chrom'].astype(int)
    df = df.sort_values(['chrom', 'pos'])

    # Calculate chromosome sizes and offsets with better spacing
    chrom_sizes = df.groupby('chrom')['pos'].max()
    max_size = chrom_sizes.max()
    spacing = max_size * 0.08  # 8% of max chromosome size for spacing
    
    # Calculate offsets with spacing
    chrom_offsets = pd.Series(0, index=range(1, 23))
    current_offset = 0
    
    for chrom in range(1, 23):
        chrom_offsets[chrom] = current_offset
        if chrom in chrom_sizes.index:
            current_offset += chrom_sizes[chrom] + spacing

    # Apply offsets to create manhattan x-coordinates
    df['x_manhattan'] = df.apply(lambda row: row['pos'] + chrom_offsets[row['chrom']], axis=1)
    
    return df, chrom_sizes, chrom_offsets

def get_y_axis_config(max_qval):
    """Determine appropriate y-axis configuration based on max -log10(p-value)."""
    if max_qval <= 14:
        tick_interval = 2
        max_range = 14
    elif max_qval <= 28:
        tick_interval = 4
        max_range = 28
    elif max_qval <= 40:
        tick_interval = 8
        max_range = 40
    elif max_qval <= 70:
        tick_interval = 10
        max_range = 70
    else:
        tick_interval = 20
        max_range = ((max_qval // 20) + 1) * 20

    return {
        'tick_interval': tick_interval,
        'max_range': max_range
    }

def create_manhattan_plot(df, chrom_sizes, chrom_offsets, output_path, height=600):
    """Create and save a manhattan plot with minimal visualization."""
    # Calculate tick positions for chromosomes
    tickvals = []
    ticktext = []
    for chrom in range(1, 23):
        if chrom in chrom_sizes.index:
            center = chrom_offsets[chrom] + chrom_sizes[chrom] / 2
            tickvals.append(center)
            ticktext.append(str(chrom))

    # Get y-axis configuration
    max_qval = df['log_p'].max()
    y_config = get_y_axis_config(max_qval)
    
    # Create plot data with alternating colors
    plot_data = []
    
    # Add significance threshold line
    # plot_data.append(go.Scatter(
    #     x=[0, df['x_manhattan'].max()],
    #     y=[-np.log10(SIGNIFICANCE_THRESHOLD), -np.log10(SIGNIFICANCE_THRESHOLD)],
    #     mode='lines',
    #     line=dict(color='gray', width=1, dash='dash'),
    #     name='Significance',
    #     showlegend=False
    # ))

    # Add chromosome data
    for chrom in range(1, 23):
        chrom_data = df[df['chrom'] == chrom]
        color_index = (chrom - 1) % 2
        
        plot_data.append(go.Scattergl(
            x=chrom_data['x_manhattan'],
            y=chrom_data['log_p'],
            mode='markers',
            marker=dict(
                color=ALTERNATING_COLORS[color_index],
                size=3,
                opacity=0.7
            ),
            hoverinfo='text',
            text=chrom_data.apply(
                lambda row: (
                    f"Chromosome: {row['chrom']}<br>"
                    f"Position: {row['pos']:,}<br>"
                    f"-log10(p): {row['log_p']:.2f}<br>"
                    f"p-value: {row['pval']:.2e}"
                ),
                axis=1
            ),
            showlegend=False
        ))

    # Create layout with minimal settings
    layout = go.Layout(
        showlegend=False,
        xaxis=dict(
            showticklabels=False,  # Hide x-axis tick labels
            showgrid=False,        # Hide grid
            zeroline=False,        # Hide zero line
            showline=False,        # Hide axis line
            ticks="",             # Hide tick marks
        ),
        yaxis=dict(
            showticklabels=False,  # Hide y-axis tick labels
            showgrid=False,        # Hide grid
            zeroline=False,        # Hide zero line
            showline=False,        # Hide axis line
            ticks="",             # Hide tick marks
        ),
        width=1500,
        height=height,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(
            l=0,    # Reduced left margin since we don't have labels
            r=0,    # Reduced right margin
            t=10,   # Reduced top margin
            b=10    # Reduced bottom margin
        ),
        hovermode='closest'
    )

    # Create and save figure
    fig = go.Figure(data=plot_data, layout=layout)
    
    # Update hover template
    fig.update_traces(
        hovertemplate="%{text}<extra></extra>"
    )

    # Save the figure
    fig.write_image(
        output_path,
        scale=3,
        width=1500,
        height=height
    )
    
    return fig
def fetch_tabix_data(chrom, gz_file_path, tbi_file_path):
    """Fetch data for a single chromosome from a tabix file."""
    results = []
    with pysam.TabixFile(gz_file_path, index=tbi_file_path) as tabix_file:
        for row in tabix_file.fetch(str(chrom)):
            fields = row.strip().split('\t')
            try:
                pos = int(float(fields[2]))
                pval = float(fields[7]) if fields[7] != 'NA' else None
                if pval is not None:
                    results.append({
                        'chrom': str(chrom),
                        'pos': pos,
                        'pval': pval,
                    })
            except (ValueError, IndexError) as e:
                print(f"Error processing row for chromosome {chrom}: {str(e)}")
                continue
    return results



def main():
    # Directory paths
    data_dir = "/nfs/platlas_stor/tabix"
    output_dir = "/nfs/platlas_stor/mh_plots"
    
    # Process each file
    for gz_file_path in glob.glob(f"{data_dir}/*pval_up_to_0.1.gz"):
        if gz_file_path.endswith('.tbi'):
            continue
            
        tbi_file_path = f"{gz_file_path}.tbi"
        if not Path(tbi_file_path).exists():
            print(f"Warning: No .tbi file found for {gz_file_path}")
            continue
            
        file_name = Path(gz_file_path).stem
        output_path = Path(output_dir) / f"manhattan_{file_name}.png"
        
        print(f"Processing {file_name}...")
        
        # Process data and create plot
        df = pd.DataFrame([
            item for chrom in range(1, 23)
            for item in fetch_tabix_data(chrom, gz_file_path, tbi_file_path)
        ])
        
        if not df.empty:
            df, chrom_sizes, chrom_offsets = process_dataframe(df)
            create_manhattan_plot(df, chrom_sizes, chrom_offsets, output_path)
            print(f"Saved plot to {output_path}")
        else:
            print(f"No valid data found in {file_name}")

if __name__ == "__main__":
    main()