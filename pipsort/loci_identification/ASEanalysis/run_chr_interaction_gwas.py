import pandas as pd
from joblib import Parallel, delayed
import statsmodels.api as sm
import sys

f_phenocovar = sys.argv[1]
f_chr_local_ancestry = sys.argv[2]

chrom = f_chr_local_ancestry.split('_')[0]
snppos = f_phenocovar.split('_')[-1]


df = pd.read_csv(f_phenocovar, sep="\t")
df2 = pd.read_csv(f_chr_local_ancestry)
df2.rename(columns={'person_id':'IID'}, inplace=True)
df3 = df.merge(df2, on='IID', how='inner')
df3.drop(columns=['FID','IID'], inplace=True)
del df2
del df

index_first_covar = df3.columns.get_loc('age')
index_last_covar = df3.columns.get_loc('snpgeno') + 1
covariates = df3.columns[index_first_covar:index_last_covar]
local_ancestry_cols = df3.columns[index_last_covar:]
covariates = list(covariates)
local_ancestry_cols = list(local_ancestry_cols)

def run_regression(col):
    # Create the interaction term
    df3['interaction'] = df3['snpgeno'] * df3[col]

    # Define the predictors (covariates + numerical column + interaction term)
    predictors = covariates + [col, 'interaction']
    X = sm.add_constant(df3[predictors])  # Add constant for the intercept
    y = df3['phenotype']

    # Fit the regression model
    model = sm.OLS(y, X)
    result = model.fit()

    # Extract the interaction term coefficient, p-value, and R-squared
    interaction_coef = result.params['interaction']
    interaction_pval = result.pvalues['interaction']
    r_squared = result.rsquared

    return {
        'gnomix_region_start': col,
        'interaction_coefficient': interaction_coef,
        'interaction_pvalue': interaction_pval,
        'r_squared': r_squared
    }

# Run the regressions in parallel
results = Parallel(n_jobs=-1)(delayed(run_regression)(col) for col in local_ancestry_cols)

# Convert results to a DataFrame for analysis
results_df = pd.DataFrame(results)

# Save the results to a file
results_df.to_csv(f'{chrom}_{snppos}_interaction_results.csv', index=False)

