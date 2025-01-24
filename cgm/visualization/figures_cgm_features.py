# plot cgm features 
import os
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
mpl.rcParams['pdf.fonttype']=42
# 
resdir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/result/cgm_meal"
os.chdir(resdir)
data=pd.read_csv('cgm_foods_manual_features.csv')
data=data[~data['foods'].isin(['Quinoa','Glucola','Glucose','WholeWheatBread'])]
plotfoods=['Beans','Berries','Bread','Grapes','Pasta','Potatoes','Rice','Rice+Fat','Rice+Fiber','Rice+Protein']
data['relative_peak']=data['peak_value']-data['baseline_glucose']
data['deltaG_at_60_mins']=data['glucose_at_60_mins']-data['baseline_glucose']
data['deltaG_at_120_mins']=data['glucose_at_120_mins']-data['baseline_glucose']
data['deltaG_at_170_mins']=data['glucose_at_170_mins']-data['baseline_glucose']
# 
features=['AUC','AUC_above_140','AUC_above_180','relative_peak','time_to_peak','AUC_above_baseline','deltaG_at_60_mins','deltaG_at_120_mins','deltaG_at_170_mins','return_to_baseline_time','slope_baseline_to_peak']
# 
def plot_bar_with_error(df, x_col, y_col, title, xlabel, ylabel,  save_fig=False, fig_name='plot.png'):
    """
    Plots a bar plot with 2*standard error bars given a DataFrame and the names of the x and y columns.
    The plot is sorted by the mean in ascending order.
    Parameters:
    df (pandas.DataFrame): the input DataFrame
    x_col (str): the name of the x column
    y_col (str): the name of the y column
    title (str): the plot title
    xlabel (str): the x-axis label
    ylabel (str): the y-axis label
    save_fig (bool): whether to save the plot or not, defaults to False
    fig_name (str): filename to save the plot as, defaults to 'plot.png'
    """
    # Calculate mean and standard error
    means = df.groupby(x_col)[y_col].mean()
    errors = df.groupby(x_col)[y_col].sem() * 2
    # Sort means by ascending order
    means = means.sort_values()
    # Plot bars and error bars
    fig, ax = plt.subplots()
    means.plot(kind='bar', yerr=errors[means.index], ax=ax, capsize=5)
    # Set axis labels and title
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # Show or save plot
    if save_fig:
        plt.savefig(fig_name,format="pdf",bbox_inches='tight')
    else:
        plt.show()
    plt.clf()
    length=len(means)
    datasv=pd.DataFrame(data={'feature': np.repeat(y_col,length),'mean':means,'se2':errors})
    if y_col=='relative_peak':
        datasv.to_csv('figext2a_'+y_col+'.csv')
        datasv2=datasv[datasv.index.isin(['Rice','Rice+Fat','Rice+Fiber','Rice+Protein'])]
        datasv2.to_csv('fig4b_'+y_col+'.csv')
    elif y_col=='AUC_above_baseline':
        datasv2=datasv[~datasv.index.isin(['Rice+Fat','Rice+Fiber','Rice+Protein'])]
        datasv2.to_csv('fig1c_'+y_col+'.csv')
    elif y_col=='time_to_peak':
        datasv2=datasv[~datasv.index.isin(['Rice+Fat','Rice+Fiber','Rice+Protein'])]
        datasv2.to_csv('fig1c_'+y_col+'.csv')
        datasv2=datasv[datasv.index.isin(['Rice','Rice+Fat','Rice+Fiber','Rice+Protein'])]
        datasv2.to_csv('fig4b_'+y_col+'.csv')
    else:
        datasv.to_csv('figext2a_'+y_col+'.csv')

data_plot=data[data['foods'].isin(plotfoods)]
plot_bar_with_error(data_plot,'foods','time_to_peak','Time_to_peak Per Food','Foods','Time to peak',save_fig=True,fig_name = 'Time_to_peak_Per_Food.pdf')
plot_bar_with_error(data_plot,'foods','relative_peak','Peak Value Per Food', 'Foods', 'Peak Value', save_fig=True, fig_name = 'Peak_Per_Food.pdf')
plot_bar_with_error(data_plot,'foods','AUC_above_baseline','AUC Per Food', 'Foods', 'AUC above baseline', save_fig=True, fig_name = 'AUC_Per_Food.pdf')
plot_bar_with_error(data_plot,'foods','return_to_baseline_time','return time Per Food','Foods','return time', save_fig=True, fig_name = 'return_time_Per_Food.pdf')

def scatter_plot(df, x_column, y_column, color_column, save_fig=False, fig_name='plot.png'):
    """
    Creates a scatter plot using the provided dataframe and columns for x and y,
    and optionally a column for color.
    Parameters:
    df (pandas.DataFrame): The dataframe to use for creating the plot.
    x_column (str): The name of the column to use for the x-axis.
    y_column (str): The name of the column to use for the y-axis.
    color_column (str): The name of the column to use for color.
    Returns:
    None
    """
    # Create the figure and axis objects
    ax = sns.scatterplot(x = x_column, y = y_column, data = df, hue = color_column)
    # Show or save plot
    if save_fig:
        plt.savefig(fig_name,format="pdf")
    else:
        plt.show()
    plt.clf()
    df2=df[[x_column,y_column,color_column]]
    df2=df2.dropna()
    df2.to_csv('figext2b_'+x_column+y_column+'.csv',index=False)

scatter_plot(data[(data['mitigator'].isna())&(data['foods'].isin(['Rice','Beans','Berries']))], 'time_to_peak', 'relative_peak','foods', save_fig=True, fig_name = 'Scatterplot_peak_and_time_to_peak.pdf')
scatter_plot(data[(data['mitigator'].isna())&(data['foods'].isin(['Rice','Beans','Berries']))], 'slope_baseline_to_peak', 'relative_peak','foods', save_fig=True, fig_name = 'Scatterplot_peak_and_slope.pdf')

def plot_food_heatmap(df,features, fig_name = 'default.png'):
    """
    Given a pandas DataFrame `df` containing columns 'food', 'feature1', 'feature2', ..., 
    group the DataFrame by 'food', compute the average of each feature for each food and normalize the features. 
    Then plot a heatmap of the resulting averages for each food and feature using seaborn.
    """
    # Group by 'food' and compute the mean of each feature
    food_averages = df.groupby('foods')[features].mean()
    # Normalize the features using StandardScaler
    scaler = StandardScaler()
    normalized_features = scaler.fit_transform(food_averages)
    normalized_df = pd.DataFrame(normalized_features, columns=food_averages.columns, index=food_averages.index)
    # Perform hierarchical clustering on rows
    dist = linkage(normalized_df.values, method='complete', metric='euclidean')
    dendrogram(dist, labels=normalized_df.index)
    plt.show()
    # Reorder rows based on clustering and plot heatmap using seaborn
    sorted_index = normalized_df.index[dendrogram(dist, no_plot=True, color_threshold=-np.inf)['leaves']]
    print(sorted_index)
    # Plot heatmap using seaborn
    sns.heatmap(normalized_df.loc[sorted_index], cmap='coolwarm')
    plt.savefig(fig_name,format="pdf",bbox_inches='tight')
    plt.clf()
    normalized_df[normalized_df.index.isin(['Beans','Berries','Bread','Grapes','Pasta','Potatoes','Rice'])].to_csv('fig1d.csv')

featuresele=['AUC','AUC_above_140','AUC_above_180','relative_peak','time_to_peak','AUC_above_baseline','deltaG_at_60_mins','deltaG_at_120_mins','deltaG_at_170_mins','return_to_baseline_time','slope_baseline_to_peak']

plot_food_heatmap(data_plot,featuresele, 'Heatmap_foods_by_features.pdf')