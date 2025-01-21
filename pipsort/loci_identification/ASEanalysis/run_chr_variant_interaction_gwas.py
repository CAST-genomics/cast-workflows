import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import statsmodels.api as sm
import sys

f_phenocovar = sys.argv[1]
f_gnomix_region_genos = sys.argv[2]
f_bim = sys.argv[3]
f_fam = sys.argv[4]

snppos = f_phenocovar.split('_')[-1]
gnomixregion = f_gnomix_region_genos.split(".")[0]


df = pd.read_csv(f_phenocovar, sep="\t")
fam = pd.read_csv(f_fam, sep="\t", header=None)
sampleids = fam[1] 
bim = pd.read_csv(f_bim, sep="\t", header=None)
varids = bim[1]
m = np.load(f_gnomix_region_genos)
m = m.T

genos = pd.DataFrame(m, columns=varids)
genos['IID'] = sampleids

df3 = df.merge(genos, on='IID', how='inner')
df3.drop(columns=['FID','IID'], inplace=True)
del df
del genos
del m

index_first_covar = df3.columns.get_loc('age')
index_last_covar = df3.columns.get_loc('snpgeno') + 1
#index_last_covar = df3.columns.get_loc('snpgeno')
covariates = df3.columns[index_first_covar:index_last_covar]
geno_cols = df3.columns[index_last_covar:]
covariates = list(covariates)
geno_cols = list(geno_cols)


def run_regression(col):
    # Create the interaction term
    df4 = df3[df3[col] >= 0].copy() #drop missing
    df4['interaction'] = df4['snpgeno'] * df4[col]

    if df4.shape[0] == 0:
        return None

    # Define the predictors (covariates + numerical column + interaction term)
    predictors = covariates + [col, 'interaction']
    #predictors = covariates + ['interaction']
    X = sm.add_constant(df4[predictors])  # Add constant for the intercept
    y = df4['phenotype']

    # Fit the regression model
    model = sm.OLS(y, X)
    result = model.fit()

    coefficients = result.params.to_dict()  # Dictionary of coefficients
    pvalues = result.pvalues.to_dict()      # Dictionary of p-values
    r_squared = result.rsquared             # R-squared value

    if col in coefficients:
        coefficients['geno'] = coefficients.pop(col)
    if col in pvalues:
        pvalues['geno'] = pvalues.pop(col)

    return {
        'varid': col,
        'coefficients': coefficients,
        'pvalues': pvalues,
        'r_squared': r_squared
    }


# Run the regressions in parallel
#results = Parallel(n_jobs=-1)(delayed(run_regression)(col) for col in local_ancestry_cols)


results = []  # Initialize an empty list to store results

# Loop through each column in local_ancestry_cols and run the regression
for col in geno_cols:
    result = run_regression(col)  # Call the function for each column
    if result:
        results.append(result)


# Convert results to a DataFrame for analysis
results_df = pd.DataFrame(results)

# Expand coefficients and p-values into separate columns
coefficients_df = pd.json_normalize(results_df['coefficients'])
pvalues_df = pd.json_normalize(results_df['pvalues'])

# Rename columns to indicate they are coefficients or p-values
coefficients_df.columns = [f'{col}_coef' for col in coefficients_df.columns]
pvalues_df.columns = [f'{col}_pval' for col in pvalues_df.columns]

# Combine expanded columns with the original results
expanded_results_df = pd.concat([results_df.drop(['coefficients', 'pvalues'], axis=1), coefficients_df, pvalues_df], axis=1)



# Save the results to a file
expanded_results_df.to_csv(f'{gnomixregion}_{snppos}_interaction_results.csv', index=False)

