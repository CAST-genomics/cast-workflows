import pandas as pd
from joblib import Parallel, delayed
import statsmodels.api as sm
import sys

f_phenocovar = sys.argv[1]
f_chr_local_ancestry = sys.argv[2]
snpposhg19 = sys.argv[3]

chrom = f_chr_local_ancestry.split('_')[0]
snppos = f_phenocovar.split('_')[-1]


df = pd.read_csv(f_phenocovar, sep="\t")
df2 = pd.read_csv(f_chr_local_ancestry)
df2.rename(columns={'person_id':'IID'}, inplace=True)

gnomix_start_regions = [int(x) for x in df2.columns[1:]]

exact_col_idx = 0
for i in range(1, len(gnomix_start_regions)):
    if int(snpposhg19) < gnomix_start_regions[i]:
        break
    exact_col_idx += 1
exact_col = str(gnomix_start_regions[exact_col_idx])
df2 = df2[['IID', exact_col]]

df3 = df.merge(df2, on='IID', how='inner')
df3.drop(columns=['FID','IID'], inplace=True)
del df2
del df
df3 = df3[df3[exact_col] != 1]
df3.loc[df3[exact_col] == 2, exact_col] = 1
df3['interaction'] = df3['snpgeno']*df3[exact_col]

index_first_covar = df3.columns.get_loc('age')
index_last_covar = df3.columns.get_loc('snpgeno')
covariates = df3.columns[index_first_covar:index_last_covar]
covariates = list(covariates)
set2 = ['snpgeno', 'interaction']


X_first = sm.add_constant(df3[covariates])
y = df3['phenotype']
model1 = sm.OLS(y, X_first).fit()
residuals = model1.resid


X_second = sm.add_constant(df3[set2])
model2 = sm.OLS(residuals, X_second).fit()

coefficients = model2.params.to_dict()  # Dictionary of coefficients
pvalues = model2.pvalues.to_dict()      # Dictionary of p-values
r_squared = model2.rsquared             # R-squared value


results = {
   'gnomix_region_start': exact_col,
   'coefficients': coefficients,
   'pvalues': pvalues,
   'r_squared': r_squared
}

# Convert results to a DataFrame for analysis
results_df = pd.DataFrame([results])

# Expand coefficients and p-values into separate columns
coefficients_df = pd.json_normalize(results_df['coefficients'])
pvalues_df = pd.json_normalize(results_df['pvalues'])

# Rename columns to indicate they are coefficients or p-values
coefficients_df.columns = [f'{col}_coef' for col in coefficients_df.columns]
pvalues_df.columns = [f'{col}_pval' for col in pvalues_df.columns]

# Combine expanded columns with the original results
expanded_results_df = pd.concat([results_df.drop(['coefficients', 'pvalues'], axis=1), coefficients_df, pvalues_df], axis=1)



# Save the results to a file
expanded_results_df.to_csv(f'{chrom}_{snppos}_interaction_results_no_het.csv', index=False)

