import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr

data = pd.read_csv('replicates.csv')
data = data.loc[~((data['subject'] == SOMESUBJECTID) & (data['dates'] == 'SOMEDATES'))]

# Define food types and variables before using them
food_types = ['Beans', 'Berries', 'Bread', 'Grapes', 'Pasta', 'Potatoes', 'Rice', 'Rice+Fat', 'Rice+Fiber', 'Rice+Protein']
variables = ['AUC_above_baseline', 'peak_value', 'baseline_glucose', 'delta_glucose']

# Initialize a dictionary to store Pearson correlation data, keyed by food type
pearson_data = {food: {} for food in food_types}



for variable in variables:
    directory = variable
    if not os.path.exists(directory):
        os.makedirs(directory)

    for food in food_types:
        food_data = data[data['foods'] == food]
        rep1 = food_data[food_data['rep'] == 1][['subject', variable]].set_index('subject')
        rep2 = food_data[food_data['rep'] == 2][['subject', variable]].set_index('subject')
        joined_data = rep1.join(rep2, lsuffix='_rep1', rsuffix='_rep2').dropna()

        if joined_data.shape[0] >= 2:
            pearson_corr = pearsonr(joined_data[f'{variable}_rep1'], joined_data[f'{variable}_rep2'])[0]

            # Store the Pearson correlation in the dictionary
            pearson_data[food][variable] = pearson_corr

            plt.figure(figsize=(8, 8))
            sns.scatterplot(x=joined_data[f'{variable}_rep1'], y=joined_data[f'{variable}_rep2'])
            plt.title(f'{food} ({variable} - Pearson Correlation: {pearson_corr:.2f})')
            plt.xlabel(f'Replicate 1 {variable}')
            plt.ylabel(f'Replicate 2 {variable}')
            all_values = pd.concat([joined_data[f'{variable}_rep1'], joined_data[f'{variable}_rep2']])
            min_val, max_val = all_values.min(), all_values.max()
            plt.xlim(min_val, max_val)
            plt.ylim(min_val, max_val)
            plt.savefig(f'{directory}/{variable}_{food}.png')
            plt.close()

# Convert the nested dictionary to a DataFrame, with food types as rows and variables as columns
pearson_df = pd.DataFrame.from_dict(pearson_data, orient='index').reset_index().rename(columns={'index': 'Food_Type'})

pearson_df = pearson_df.round(2)

# Save the DataFrame to a CSV file
pearson_df.to_csv('pearson_coefficients.csv', index=False)
