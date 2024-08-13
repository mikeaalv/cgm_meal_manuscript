import pandas as pd
from Regression_analysis.data.variables import *

def load_and_preprocess_cgm_foods():
    """
    Load and preprocess the cgm_foods data.
    """
    cgm_foods = pd.read_csv('Regression_analysis/data/cgm_foods_manual_features.csv')
    cgm_foods = cgm_foods[(cgm_foods['food'] != 'Glucose') & (cgm_foods['food'] != 'Quinoa')]
    print(cgm_foods.columns)

    cgm_foods = pd.get_dummies(cgm_foods, columns=['food', 'mitigator'])
    cgm_foods = cgm_foods[['subject', 'AUC_above_baseline', 'AUC_above_baseline_120mins','peak_value','baseline_glucose','rep','dates'] + [col for col in cgm_foods.columns if 'food' in col or 'mitigator' in col]]
    cgm_foods['delta_glucose'] = cgm_foods['peak_value'] - cgm_foods['baseline_glucose']

    # Extract the columns that correspond to food types
    food_mitigator_columns = [col for col in cgm_foods.columns if 'food' in col or 'mitigator' in col]
    
    # Group by subject and food type and compute the mean for the desired columns (delete this if we want for replicates)
    grouped = cgm_foods.groupby(['subject'] + food_mitigator_columns).agg({
        'AUC_above_baseline': 'mean',
        'AUC_above_baseline_120mins': 'mean',
        'peak_value': 'mean',
        'baseline_glucose': 'mean',
        'delta_glucose': 'mean'
    }).reset_index()

    # Group by subject and select the first occurrence (for predictive_clinical.csv)
    grouped = cgm_foods.groupby('subject').first().reset_index()

    # Only one food
    #grouped = grouped[grouped['food_Potatoes'] & ~grouped['mitigator_Fiber'] & ~grouped['mitigator_Protein'] & ~grouped['mitigator_Fat']]

    return grouped

def load_and_metabolomics():
    # Load the data
    metabolomics_data = pd.read_csv('Regression_analysis/data/metabolomics.csv')
    column_names = metabolomics_data.columns
    parsed_columns = [name.split("_") for name in column_names]

    alphanumeric_parts = [item[0] for item in parsed_columns]
    date_parts = [item[1] for item in parsed_columns]

    numeric_parts = [int(part[1:]) for part in alphanumeric_parts if part[1:].isdigit()]

    date_objects = pd.to_datetime(date_parts, errors='coerce')

    metabolomics_data = metabolomics_data.T
    metabolomics_data['date'] = date_objects
    metabolomics_data['subject'] = numeric_parts

    metabolomics_data = metabolomics_data.groupby('subject').mean().reset_index()

    cols =  [col for col in metabolomics_data.columns if col not in ['subject', 'date']] + ['subject', 'date']
    metabolomics_data = metabolomics_data[cols]

    columns_copy = metabolomics_data.columns.tolist()
    columns_copy[:-2] = ['meta_' + str(col) for col in columns_copy[:-2]]
    metabolomics_data.columns = columns_copy

    metabolomics_data.to_csv('Regression_analysis/data/metabolomics_processed.csv', index=False)
    return metabolomics_data

def load_and_lipidomics():
    # Load the data
    lipidomics_data = pd.read_csv('Regression_analysis/data/lipidomics.csv')
    column_names = lipidomics_data.columns
    parsed_columns = [name.split("_") for name in column_names]

    alphanumeric_parts = [item[0] for item in parsed_columns]
    date_parts = [item[1] for item in parsed_columns]

    numeric_parts = [int(part[1:]) for part in alphanumeric_parts if part[1:].isdigit()]

    date_objects = pd.to_datetime(date_parts, errors='coerce')

    lipidomics_data = lipidomics_data.T
    lipidomics_data['date'] = date_objects
    lipidomics_data['subject'] = numeric_parts

    lipidomics_data = lipidomics_data.groupby('subject').mean().reset_index()

    cols =  [col for col in lipidomics_data.columns if col not in ['subject', 'date']] + ['subject', 'date']
    lipidomics_data = lipidomics_data[cols]

    columns_copy = lipidomics_data.columns.tolist()
    columns_copy[:-2] = ['lip_' + str(col) for col in columns_copy[:-2]]
    lipidomics_data.columns = columns_copy

    lipidomics_data.to_csv('Regression_analysis/data/lipidomics_processed.csv', index=False)
    return lipidomics_data

def load_and_proteomics():
    # Load the data
    proteomics_data = pd.read_csv('Regression_analysis/data/proteomics.csv')
    column_names = proteomics_data.columns

    alphanumeric_parts = column_names

    numeric_parts = [int(part[1:]) for part in alphanumeric_parts if part[1:].isdigit()]


    proteomics_data = proteomics_data.T
    proteomics_data['subject'] = numeric_parts
    
    
    cols =  [col for col in proteomics_data.columns if col not in ['subject']] + ['subject']
    proteomics_data = proteomics_data[cols]

    columns_copy = proteomics_data.columns.tolist()
    columns_copy[:-1] = ['prot_' + str(col) for col in columns_copy[:-1]]
    proteomics_data.columns = columns_copy

    proteomics_data.to_csv('Regression_analysis/data/proteomics_processed.csv', index=False)
    return proteomics_data

def merge_data(cgm_foods, metabolomics, lipidomics, proteomics, metadata):
    """
    Merge the cgm_foods and metabolomics dataframes and save to CSV.
    """
    print(f"CGM: {cgm_foods['subject'].nunique()}")
    df = pd.merge(cgm_foods, metabolomics, on='subject', how='outer')
    print(f"metabolomics: {df['subject'].nunique()}")
    df = pd.merge(df, lipidomics, on='subject',how='outer')
    print(f"lip: {df['subject'].nunique()}")
    df = pd.merge(df, proteomics, on='subject', how='outer')
    print(f"prot: {df['subject'].nunique()}")
    df = pd.merge(df, metadata,  on = 'subject',how='outer')
    print(f"metabolomics: {df['subject'].nunique()}")

    return df

def load_and_preprocess_metadata(columns_to_keep):
    """
    Load and preprocess the metadata.
    """
    metadata = pd.read_csv('Regression_analysis/data/metadata_clean_all_ver2.csv')
    metadata['sex'] = metadata['sex.factor']
    metadata['subject'] = metadata['study_id'].str[-3:].astype(int)
    metadata = metadata.drop(columns=['study_id'])
    #metadata = metadata.dropna(subset = ['sex'])
    metadata['sex'] = metadata['sex'].replace({'Male': 1, 'Female': 0})
    metadata['ethnicity'] = metadata['ethnicity'].replace({'Asian, White':'Asian','Hispanic or Latino':'HispanicLatino'})
    metadata = pd.get_dummies(metadata, columns=['ethnicity'])
    metadata = metadata[columns_to_keep+['subject']]
    return metadata



def main():
    """
    Main function to preprocess data
    """
    cgm_foods = load_and_preprocess_cgm_foods()
    metabolomics = load_and_metabolomics()
    lipidomics = load_and_lipidomics()
    proteomics = load_and_proteomics()
    metadata = load_and_preprocess_metadata(columns_to_keep=['a1c_avg_all',"sspg_avg_all","ogtt_t_120_avg_all","bmi_avg_all", "ie_heyjun", "modified_DI_heyjun", "insulin_fasting_avg_all","fbg_avg_all", "hepatic_IR_heyjun", "ffa_avg_heyjun", "age_today_avg_all", "sex"])

    df = merge_data(cgm_foods, metabolomics, lipidomics, proteomics, metadata)
    print(cgm_foods.head())
    print(metabolomics.head())
    print(metadata.head())
    #df.to_csv('Regression_analysis/data/predictive_clinical_ver2.csv')
    df.to_csv('Regression_analysis/data/predictive_data_ver4.csv')
    # Can create replicates.csv, predictive_data_final.csv, or predictive_clinical.csv
    # At the end for predictive_clinical.  There is 39 participants which seems about right



if __name__ == "__main__":
    main()