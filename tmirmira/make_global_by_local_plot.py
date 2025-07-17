import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys

chrom = sys.argv[1]
local_a = sys.argv[2]

ancestrydf = pd.read_csv("ancestry_info.csv")
ancestrydf = ancestrydf.rename(columns={'research_id':'person_id'})
codes = {'eur':0,'eas':1,'nat':2,'afr':3,'sas':4,'mid':7}

#color map
base_colors = {
    'afr': '#d62728',  # red
    'eur': '#1f77b4',  # blue
    'amr': '#2ca02c',  # green
    'eas': '#ff7f0e',  # orange
    'mid': '#9467bd',  # purple
    'sas': '#8c564b'   # brown
}

# Function to lighten a color
def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    r, g, b = mc.to_rgb(c)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    return mc.to_hex(colorsys.hls_to_rgb(h, min(1, l + amount * (1 - l)), s))

# Build full palette
ancestries = ['afr', 'eur', 'amr', 'eas', 'mid', 'sas']
color_map = {}

for anc in ancestries:
    color_map[anc] = base_colors[anc]
    color_map[f'{anc}_oth'] = lighten_color(base_colors[anc], amount=0.4)


def plot_percent(chrom, local_a):
    gnomix = f'gnomix-chr{chrom}.msp'
    m = np.loadtxt(gnomix, skiprows=2)

    with open(gnomix, 'r') as file:
        second_line = file.readlines()[1]
    array = np.array(second_line.strip().split('\t'))
    columns = list(array)
    df = pd.DataFrame(m, columns=columns)
    del m
    columns = df['spos'].astype('int')
    df = df.T
    df = df[6:]
    df.columns = columns

    local_a_code = codes[local_a]
    df = df==local_a_code
    df = df.astype('int')
    df_pair_sum = df.iloc[::2].copy()
    df_pair_sum.iloc[:, :] += df.iloc[1::2].values

    # Extract person_id from the original index
    df_pair_sum['person_id'] = df.index[::2].map(lambda x: str(x).split('.')[0])
    df_pair_sum = df_pair_sum.reset_index(drop=True)
    df_pair_sum['person_id'] = df_pair_sum['person_id'].astype('int64')
    local_and_global = df_pair_sum.merge(ancestrydf[['person_id', 'ancestry_my_oth']], on='person_id', how='inner')
    local_and_global = local_and_global.drop(columns='person_id')

    pos_cols = [col for col in local_and_global.columns if col != 'ancestry_my_oth']

    plt.figure(figsize=(11, 6))

    # Step 3: For each ancestry, compute frequencies
    for anc in local_and_global['ancestry_my_oth'].unique():
        sub = local_and_global[local_and_global['ancestry_my_oth'] == anc]
        denom = 2 * len(sub)  # max value if all were '2'
        perc = sub[pos_cols].sum(axis=0) / denom
        plt.plot(pos_cols, perc, label=anc, color=color_map.get(anc, 'gray'))

    # Step 5: Plot
    plt.title(f"Chromosome {chrom}")
    plt.xlabel("Genomic Position")
    plt.ylim(-0.05,1.05)
    plt.ylabel(f"% {local_a} ancestry")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"chrom{chrom}_{local_a}_local_ancestry.png")
    plt.show()

plot_percent(chrom, local_a)
