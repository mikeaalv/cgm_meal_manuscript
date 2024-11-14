import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42

#data = pd.read_csv('replicates2.csv')
data = pd.read_csv('../../data/cgm_foods_manual_features.csv')
data = data.loc[~((data['subject'] == "SOMEID") & (data['dates'] == 'SOMEDATE'))]

food_types = ['Beans', 'Berries', 'Bread', 'Grapes', 'Pasta', 'Potatoes', 'Rice', 'Rice+Fat', 'Rice+Fiber', 'Rice+Protein']
variables = ['AUC_above_baseline', 'AUC_above_baseline_120mins', 'peak_value', 'baseline_glucose']

pearson_data = {food: {} for food in food_types}

for variable in variables:
    directory = variable
    if not os.path.exists(directory):
        os.makedirs(directory)
    fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 8))
    axs = axs.flatten()

    global_min = data[variable].min()
    global_max = data[variable].max()

    # Add a buffer around the edges
    buffer_size = 0.05 * (global_max - global_min)  # 5% buffer
    buffered_min = global_min - buffer_size
    buffered_max = global_max + buffer_size

    for idx, food in enumerate(food_types):
        food_data = data[data['foods'] == food]
        rep1 = food_data[food_data['rep'] == 1][['subject', variable]].set_index('subject')
        rep2 = food_data[food_data['rep'] == 2][['subject', variable]].set_index('subject')
        joined_data = rep1.join(rep2, lsuffix='_rep1', rsuffix='_rep2').dropna()

        if joined_data.shape[0] >= 2:
            pearson_corr = pearsonr(joined_data[f'{variable}_rep1'], joined_data[f'{variable}_rep2'])[0]
            pearson_data[food][variable] = pearson_corr

            ax = axs[idx]
            sns.scatterplot(x=joined_data[f'{variable}_rep1'], y=joined_data[f'{variable}_rep2'], ax=ax)
            ax.set_title(f'{food} (r={pearson_corr:.2f})')
            ax.set_xlabel('Replicate 1')
            ax.set_ylabel('Replicate 2')

            # Set axis limits with buffer
            ax.set_xlim(buffered_min, buffered_max)
            ax.set_ylim(buffered_min, buffered_max)

    plt.tight_layout()
    plt.savefig(f'{directory}/{variable}_all_foods.pdf', format = 'pdf', bbox_inches = 'tight')
    plt.close()

pearson_df = pd.DataFrame.from_dict(pearson_data, orient='index').reset_index().rename(columns={'index': 'Food_Type'})
pearson_df = pearson_df.round(2)
pearson_df.to_csv('pearson_coefficients.csv', index=False)





variables = ['AUC_above_baseline', 'AUC_above_baseline_120mins']

# Initialize dictionary to store Pearson correlations
pearson_data = {food: {var: np.nan for var in variables} for food in food_types+['all_standardized_foods']}

all_standardized_foods_AUC = []#{var:[] for var in variables}
all_standardized_foods_AUC_120mins = []
# Calculate Pearson correlations
for food in food_types:
    for variable in variables:
        food_data = data[data['foods'] == food]
        rep1 = food_data[food_data['rep'] == 1][['subject', variable]].set_index('subject')
        rep2 = food_data[food_data['rep'] == 2][['subject', variable]].set_index('subject')
        joined_data = rep1.join(rep2, lsuffix='_rep1', rsuffix='_rep2').dropna()

        if len(joined_data) >= 2:
            pearson_corr = pearsonr(joined_data[variable + '_rep1'], joined_data[variable + '_rep2'])[0]
            pearson_data[food][variable] = pearson_corr
        if food in ['Beans', 'Berries', 'Bread', 'Grapes', 'Pasta', 'Potatoes', 'Rice']:
            if variable == 'AUC_above_baseline':
                all_standardized_foods_AUC.append(joined_data)
            else:
                all_standardized_foods_AUC_120mins.append(joined_data)

all_standardized_foods_AUC_data = pd.concat(all_standardized_foods_AUC, ignore_index=True)
pearson_corr = pearsonr(all_standardized_foods_AUC_data['AUC_above_baseline_rep1'], all_standardized_foods_AUC_data['AUC_above_baseline_rep2'])[0]
pearson_data['all_standardized_foods']['AUC_above_baseline'] = pearson_corr

all_standardized_foods_AUC_120mins_data = pd.concat(all_standardized_foods_AUC_120mins, ignore_index=True)
pearson_corr = pearsonr(all_standardized_foods_AUC_120mins_data['AUC_above_baseline_120mins_rep1'], all_standardized_foods_AUC_120mins_data['AUC_above_baseline_120mins_rep2'])[0]
pearson_data['all_standardized_foods']['AUC_above_baseline_120mins'] = pearson_corr    

# Create DataFrame from dictionary
pearson_df = pd.DataFrame.from_dict(pearson_data, orient='index')

# Update column names and label
new_labels = {'AUC_above_baseline': 'AUC_above_baseline_180mins', 'AUC_above_baseline_120mins': 'AUC_above_baseline_120mins'}
pearson_df.columns = [new_labels.get(var, var) for var in pearson_df.columns]

pearson_df.to_csv('replicates_pearson.csv')

# Plotting
fig, ax = plt.subplots(figsize=(14, 6))
pearson_df.plot(kind='bar', ax=ax, rot=0)
ax.set_title('Pearson Correlations for AUC Metrics Across Food Types')
ax.set_ylabel('Pearson Correlation')
ax.set_xlabel('Food Type')
plt.xticks(rotation=45)
plt.legend(title="Metrics")
plt.tight_layout()
plt.savefig('replicates_bar_chart.pdf', format = 'pdf', bbox_inches = 'tight')
