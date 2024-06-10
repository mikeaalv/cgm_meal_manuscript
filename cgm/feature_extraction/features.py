#features.py From Ben

import pandas as pd
import numpy as np
from scipy.integrate import trapezoid


# Return to baseline from peak and from baseline
# Slope from peak back to baseline
# Secondary bump in glucose? # Multiple peaks?  True/False. 


# Helper functions
def post_event_filter(cgm_df, glucose_col = 'glucose'):
    return cgm_df[cgm_df['mins_since_start']>=0]

# Generic glucose features
def AUC(cgm_df, glucose_col = 'glucose'):
    return trapezoid(cgm_df[glucose_col],dx = 5)

def AUC_above_x(cgm_df, x, glucose_col = 'glucose'):
    new_df = cgm_df.copy()
    new_df[glucose_col] = cgm_df[glucose_col]-x
    new_df.loc[new_df[glucose_col]<0,glucose_col] = 0
    return AUC(new_df)

def AUC_above_140(cgm_df, glucose_col = 'glucose'):
    return AUC_above_x(cgm_df,140,glucose_col=glucose_col)

def AUC_above_180(cgm_df, glucose_col = 'glucose'):
    return AUC_above_x(cgm_df,180,glucose_col=glucose_col)

def CV(cgm_df, glucose_col = 'glucose'):
    return np.std(cgm_df[glucose_col]) / np.mean(cgm_df[glucose_col])

def average_glucose(cgm_df, glucose_col = 'glucose'):
    return np.mean(cgm_df[glucose_col])

# Excursion specific features

def AUC_after_event(cgm_df, glucose_col = 'glucose'):
    return AUC(cgm_df)

def total_length(cgm_df, glucose_col = 'glucose'):
    return len(cgm_df)

def length_after_0min(cgm_df, glucose_col = 'glucose'):
    return len(cgm_df[cgm_df['mins_since_start']>=0])

def peak_value(cgm_df, glucose_col = 'glucose'):
    return np.max(cgm_df[glucose_col])

def time_to_peak(cgm_df, glucose_col = 'glucose'):
    peak_idx = np.argmax(cgm_df[glucose_col])
    return cgm_df.iloc[peak_idx]['mins_since_start']

def baseline_glucose(cgm_df, glucose_col = 'glucose'):
    idx = (abs(cgm_df['mins_since_start'])).idxmin()
    return cgm_df.loc[idx][glucose_col]

def AUC_above_baseline(cgm_df, glucose_col = 'glucose'):
    return AUC_above_x(cgm_df,baseline_glucose(cgm_df,glucose_col=glucose_col),glucose_col=glucose_col)

def glucose_at_x_mins(cgm_df, x, glucose_col = 'glucose'):
    idx = (abs(cgm_df['mins_since_start']-x)).idxmin()
    return cgm_df.loc[idx][glucose_col]

def glucose_at_60_mins(cgm_df, glucose_col = 'glucose'):
    return glucose_at_x_mins(cgm_df, 60, glucose_col=glucose_col)

def glucose_at_120_mins(cgm_df, glucose_col = 'glucose'):
    return glucose_at_x_mins(cgm_df, 120, glucose_col=glucose_col)

def glucose_at_170_mins(cgm_df, glucose_col = 'glucose'):
    return glucose_at_x_mins(cgm_df, 170, glucose_col=glucose_col)

def return_to_baseline_time(cgm_df, glucose_col = 'glucose'):
    bg = baseline_glucose(cgm_df,glucose_col=glucose_col)
    peak_time = time_to_peak(cgm_df,glucose_col=glucose_col)
    new_df = cgm_df[cgm_df['mins_since_start']>peak_time]
    below_baseline = new_df[new_df[glucose_col]<=bg]
    if len(below_baseline)==0:
        return 170
    else:
        return below_baseline.iloc[0]['mins_since_start']

def return_to_baseline(cgm_df, glucose_col = 'glucose'):
    return return_to_baseline_time(cgm_df)<170

def slope_baseline_to_peak(cgm_df, glucose_col = 'glucose'):
    bg = baseline_glucose(cgm_df, glucose_col=glucose_col)
    peak_glucose = peak_value(cgm_df, glucose_col=glucose_col)
    denominator = time_to_peak(cgm_df, glucose_col=glucose_col)
    return (peak_glucose-bg)/denominator

# Call all functions
def calculate_all_features(cgm_df,glucose_col = 'glucose'):
    functions = {
        'whole_time_to_peak':time_to_peak,
        'total_length':total_length,
        'length_after_0min':length_after_0min
        }
    return {f:functions[f](cgm_df) for f in functions.keys()}

def calculate_all_after_event_features(cgm_df, glucose_col = 'glucose'):
    functions = {
        'AUC':AUC,
        'AUC_above_140':AUC_above_140,
        'AUC_above_180':AUC_above_180,
        'CV':CV,
        'average_glucose':average_glucose,
        'peak_value':peak_value,
        'time_to_peak':time_to_peak,
        'baseline_glucose':baseline_glucose,
        'AUC_above_baseline':AUC_above_baseline,
        'glucose_at_60_mins':glucose_at_60_mins,
        'glucose_at_120_mins':glucose_at_120_mins,
        'glucose_at_170_mins':glucose_at_170_mins,
        'return_to_baseline_time':return_to_baseline_time,
        'return_to_baseline':return_to_baseline,
        'slope_baseline_to_peak':slope_baseline_to_peak
    }
    return {f:functions[f](cgm_df) for f in functions.keys()}

def calculate_food_challenge_features(cgm_df,glucose_col = 'glucose'):
    result = cgm_df.groupby(['subject','foods','rep','food','mitigator'], dropna=False).apply(calculate_all_features)
    result = result.apply(pd.Series).reset_index()
    filtered_df = post_event_filter(cgm_df, glucose_col=glucose_col)
    post_event_result = filtered_df.groupby(['subject','foods','rep','food','mitigator'], dropna=False).apply(calculate_all_after_event_features)
    post_event_result = post_event_result.apply(pd.Series).reset_index()
    return pd.merge(result,post_event_result, on =['subject','foods','rep','food','mitigator'])

