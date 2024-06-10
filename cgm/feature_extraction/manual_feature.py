import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
# 
import importlib.util
import sys
import os
# spec=importlib.util.spec_from_file_location("features", "/Users/yuewu/Documents/GitHub/cgm_analysis/scripts/features.py")
# features=importlib.util.module_from_spec(spec)
# sys.modules["features"]=features
# spec.loader.exec_module(features)
import features
# 
respath="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/result/cgm_meal/"
os.chdir(respath)
cgm_foods=pd.read_csv(respath+"cgm_foods_smooth.txt",sep='\s+')
# 
food_summary_statistics=features.calculate_food_challenge_features(cgm_foods)
food_summary_statistics.to_csv('cgm_foods_manual_features.csv',index=False)
# 
food_summary_statistics2=food_summary_statistics[food_summary_statistics['foods'].isin(['Rice','Beans'])]
groups=food_summary_statistics2.groupby('food')
fig,ax=plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name,group in groups:
    ax.plot(group.slope_baseline_to_peak,group.peak_value,marker='o',linestyle='',ms=2,label=name)
ax.legend()
plt.show()
