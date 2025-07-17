import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
import sys

chr1 = sys.argv[1]
chr2 = sys.argv[2]
local_a1 = sys.argv[3]
local_a2 = sys.argv[4]

codes = {'eur':0,'eas':1,'nat':2,'afr':3,'sas':4,'mid':7}

def corr_mat_v3(chr1, chr2, local_a1, local_a2):

    #read chr1
    gnomix1 = f'gnomix-chr{chr1}.msp'
    m = np.loadtxt(gnomix1, skiprows=2)

    with open(gnomix1, 'r') as file:
        second_line = file.readlines()[1]
    array = np.array(second_line.strip().split('\t'))
    columns = list(array)
    df = pd.DataFrame(m, columns=columns)
    del m
    columns = df['spos'].astype('int')
    df = df.T
    df = df[6:]
    df.columns = columns

    local_a_code = codes[local_a1]
    df = df==local_a_code
    df1 = df.astype('int')
    df_pair_sum1 = df.iloc[::2].copy()
    df_pair_sum1.iloc[:, :] += df.iloc[1::2].values
    df_pair_sum1 = df_pair_sum1>0

    #read chr2
    gnomix2 = f'gnomix-chr{chr2}.msp'
    m = np.loadtxt(gnomix2, skiprows=2)

    with open(gnomix2, 'r') as file:
        second_line = file.readlines()[1]
    array = np.array(second_line.strip().split('\t'))
    columns = list(array)
    df = pd.DataFrame(m, columns=columns)
    del m
    columns = df['spos'].astype('int')
    df = df.T
    df = df[6:]
    df.columns = columns

    local_a_code = codes[local_a2]
    df = df==local_a_code
    df2 = df.astype('int')
    df_pair_sum2 = df.iloc[::2].copy()
    df_pair_sum2.iloc[:, :] += df.iloc[1::2].values
    df_pair_sum2 = df_pair_sum2>0



    X = df_pair_sum1.to_numpy()
    Y = df_pair_sum2.to_numpy()

    # Center the columns (subtract mean)
    X_centered = X - X.mean(axis=0)
    Y_centered = Y - Y.mean(axis=0)

    # Compute standard deviations
    X_std = X_centered.std(axis=0)
    Y_std = Y_centered.std(axis=0)

    # Avoid division by zero by replacing 0 std with 1 temporarily
    X_std_safe = np.where(X_std == 0, 1, X_std)
    Y_std_safe = np.where(Y_std == 0, 1, Y_std)

    # Normalize
    X_norm = X_centered / X_std_safe
    Y_norm = Y_centered / Y_std_safe

    # Set columns with std = 0 to all zeros
    X_norm[:, X_std == 0] = 0
    Y_norm[:, Y_std == 0] = 0

    correlation_matrix = X_norm.T @ Y_norm / X.shape[0]  # or np.dot(X_norm.T, Y_norm) / X.shape[0]

    correlation_matrix = correlation_matrix.astype(np.float32)
    np.save(f"chrom{chr1}_chrom{chr2}_{local_a1}_vs_{local_a2}.npy", correlation_matrix)

    plt.figure(figsize=(16, 9))
    sns.heatmap(correlation_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title(f"{local_a1} vs {local_a2} local ancestry")
    plt.xlabel(f"Chrom {chr2}")
    plt.ylabel(f"Chrom {chr1}")
    plt.tight_layout()
    plt.savefig(f"chrom{chr1}_chrom{chr2}_{local_a1}_vs_{local_a2}.png")
    plt.show()

corr_mat_v3(chr1,chr2,local_a1,local_a2)
