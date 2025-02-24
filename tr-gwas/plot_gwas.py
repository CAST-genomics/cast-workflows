"""
Make Manhattan and QQ plots
QQ plot code copied from:
https://github.com/satchellhong/qqman/
"""

import matplotlib.pyplot as plt
import numbers
import numpy as np
import pandas as pd
import seaborn as sns
import argparse
import os

#download summary stats file from workspace bucket
# gsutil cp {WORKSPACE_BUCKET}/tr-gwas_result/

def GetOutPath(file_path):
    return file_path.replace(".tab", "")

def read_file(file_path):  
    data = []
    # Initialize a flag for whether we are reading data after the header
    is_reading_data = False
    columns = []
    header_found = False  # To track if we've already found the first header
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    # Open the file
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespace
            if line.startswith('==>'):
                continue
            if line.startswith('#CHROM'):
                if not header_found:
                    columns = line.split()  # Store the first header
                    header_found = True  # Mark that we've found the first header
                continue  # Skip all subsequent header lines

            if header_found:
                # Only process non-empty lines
                if line:
                    data.append(line.split())  # Split the line by whitespace or tab
    df = pd.DataFrame(data, columns=columns)
    return df


def PlotManhattan(df, outpath):
    #drop na values
    df = df[df['TEST'] == 'ADD']
    df['P'] = df['P'].replace('NA', np.nan)
    df_cleaned = df.dropna(subset=['P'])
    print(f'{df_cleaned.shape[0]} number of STRs are left after removing NA')
    #convert str to numerical value 
    df_cleaned['P'] = df_cleaned['P'].astype(float)
    df_cleaned['#CHROM'] = df_cleaned['#CHROM'].astype(int)
    df_cleaned['POS'] = df_cleaned['POS'].astype(int)
    df_cleaned['-log10p'] = -np.log10(df_cleaned.P)
    df_sorted = df_cleaned.sort_values(by=["#CHROM", "POS"])
    df_sorted.reset_index(inplace=True, drop=True)
    # save clean sorted dataframe
    df_sorted.save_csv(f'{outpath}_ADD.csv',index=False)
    df_sorted['i'] = df_sorted.index
    # Generate Manhattan plot
    plot = sns.relplot(data=df_sorted, x='i', y='-log10p', s=6, aspect=4, linewidth=0,
                       hue='#CHROM', palette="tab10", legend=None)
    chrom_df = df_sorted.groupby('#CHROM')['i'].median()
    plot.ax.set_xlabel('Chromosomes')
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle(f'Manhattan plot of STR genome-wide association on AllofUs {outpath}')
    plot.ax.axhline(8, linestyle='--', linewidth=1)
    plot.savefig(f'{outpath}.manhattan.png', dpi=700)

def ppoints(n, a=None):
    """ numpy analogue of `R`'s `ppoints` function """
    if isinstance(n, numbers.Number):
        n = float(n)
    else:
        n = float(len(n))
    if a is None:
        a = .375 if n <= 10 else .5
    return (np.arange(n) + 1 - a) / (n + 1 - 2 * a)

def PlotQQ(df, outpath):
    df['P'] = df['P'].replace('NA', np.nan)
    df = df.dropna(subset=['P'])
    df['P'] = df['P'].astype(float)
    df['-log10p'] = -np.log10(df.P)
    p_vals = df["-log10p"].dropna().apply(lambda x: 10**-x)
    p_vals = p_vals[(0 < p_vals) & (p_vals < 1)]
    observed = -np.log10(np.sort(np.array(p_vals)))
    expected = -np.log10(ppoints(len(p_vals)))
    x_padding = (np.nanmax(expected) - np.nanmin(expected)) / 12
    y_padding = (np.nanmax(observed) - np.nanmin(observed)) / 12

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 12))

    ax.scatter(expected, observed, c='k')

    xlim_min = np.nanmin(expected) - x_padding
    xlim_max = np.nanmax(expected) + x_padding
    ylim_min = np.nanmin(observed) - y_padding
    ylim_max = np.nanmax(observed) + y_padding

    max_lim = xlim_max if xlim_max < ylim_max else ylim_max
    min_lim = xlim_min if xlim_min > ylim_min else ylim_min
    ax.plot([min_lim, max_lim], [min_lim, max_lim], 'r-')

    ax.set_xlim([xlim_min, xlim_max])
    ax.set_ylim([ylim_min, ylim_max])
    ax.set_xlabel("Expected $-log_{10}(p)$")
    ax.set_ylabel("Observed $-log_{10}(p)$")

    fig.savefig(outpath+".qq.png")

def main():
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots for GWAS data.")
    parser.add_argument("--gwas", help="GWAS summary stats file path", type=str, required=True)
    args = parser.parse_args()


    # Create output directory if it does not exist
    outdir = "outputs"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    
    # Read data
    gwas_df = read_file(args.gwas)  
    filename = GetOutPath(args.gwas)
    outpath = f'{outdir}/{filename}'
    PlotManhattan(gwas_df, outpath)
    PlotQQ(gwas_df, outpath)

if __name__ == "__main__":
    main()
