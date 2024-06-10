# initial preprocessing of the cgm data
import os
import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import timedelta
# path
respath="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/result/cgm_meal_preprocessing/"
cgmdatapath="/PATHTOAWDATA/"
proc_para_path="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/data/cgm_preprocessing_information/"
# food_path=rawdatapath+"cronometer/"
os.chdir(respath)
# 
def read_new_cgm(cgmdatapath,subject,file):
    new_cgm=pd.read_csv(cgmdatapath+file)
    new_cgm['datetime_local']=new_cgm['Timestamp (YYYY-MM-DDThh:mm:ss)']
    new_cgm['glucose']=new_cgm['Glucose Value (mg/dL)']
    new_cgm['subject']=subject
    new_cgm=new_cgm[['datetime_local','glucose','subject']]
    new_cgm=new_cgm.dropna(subset='datetime_local')
    new_cgm['datetime_local']=pd.to_datetime(new_cgm['datetime_local'])
    return new_cgm
def read_new_cgm_txt(cgmdatapath,subject,file):
    new_cgm=pd.read_table(cgmdatapath+file,index_col=False)
    new_cgm['datetime_local']=new_cgm['GlucoseDisplayTime']
    new_cgm['glucose']=new_cgm['GlucoseValue']
    new_cgm['subject']=subject
    new_cgm=new_cgm[['datetime_local','glucose','subject']]
    new_cgm=new_cgm.dropna(subset='datetime_local')
    new_cgm['datetime_local']=pd.to_datetime(new_cgm['datetime_local'])
    # new_cgm['calibration']=False
    # 
    # new_cgm2=pd.DataFrame()
    # new_cgm2['datetime_local']=new_cgm['MeterDisplayTime']
    # new_cgm2['glucose']=new_cgm['MeterValue']
    # new_cgm2['subject']=subject
    # new_cgm2['calibration']=True
    # new_cgm=pd.concat([new_cgm,new_cgm2],ignore_index=True)
    return new_cgm

# metadata 
cgm_df=pd.read_csv(proc_para_path+'all_cgm_data_dates.csv',parse_dates=['time'])
cgm_df['cgm']=pd.to_numeric(cgm_df['cgm'],errors='coerce')
cgm_df['subject']=pd.to_numeric(cgm_df['id'].str[-3:])
cgm_df=cgm_df[['time','cgm','subject','cycle']].sort_values(['subject','time'])
cgm_df=cgm_df[pd.notnull(cgm_df['time'])]
cgm_df=cgm_df.rename({
    'time':'datetime_local',
    'cgm':'glucose'
},axis=1)
# new data
tab_newdata=pd.read_csv(proc_para_path+"newadd_files.csv")
new_data=pd.DataFrame()
for ind,row in tab_newdata.iterrows():
    if row[0]=="txt":
        loctab=read_new_cgm_txt(cgmdatapath,row[1],row[2])
    else:
        loctab=read_new_cgm(cgmdatapath,row[1],row[2])
    new_data=pd.concat([new_data, loctab])
## 
########### SUBJECT specific loading that's masked#################
# 
data=pd.concat([cgm_df,new_data],ignore_index=True)
data=data[['datetime_local','glucose','subject']]
data['glucose']=pd.to_numeric(data['glucose'],errors='coerce')
data=data.drop_duplicates()
data=data.sort_values(['subject','datetime_local'])
data.to_csv('final_cgm.csv',index=False)
# Extract food log info
food_log=pd.DataFrame()
cycles_info=pd.DataFrame()
for  root,dirs,files in os.walk(cgmdatapath):
    for fileele in files:
        cron_mask=".cron" in fileele.lower()
        cron_otherinfr="exercise" in fileele.lower() or "biometr" in fileele.lower() or "ogtt" in fileele.lower()
        if cron_mask and not cron_otherinfr:
            
            id=root.split("/")[7].split("_")[0]
            cycle=root.split("/")[-1]
            try:
                data=pd.read_csv(root+"/"+fileele,parse_dates=[['Day','Time']])
            except:
                print(root+"/"+fileele)
                print('continued! couldnt read')
                continue
            if 'Unnamed' in data.columns[1]:
                data.rename(columns={data.columns[1]:'Group'}, inplace=True)
            if 'Group' not in data.columns or 'Food Name' not in data.columns:
                print(root+"/"+fileele)
                print('continued! no Group')
                continue
            data=data[['Day_Time','Group','Food Name']]
            data['id']=id
            data['cycle']=cycle
            food_log=pd.concat([food_log,data],ignore_index=True)
            break
food_log.to_csv('new_food.csv')
# Basic cleaning
food_log['Day_Time']=pd.to_datetime(food_log['Day_Time'],errors='coerce')
food_log=food_log.dropna(subset=['Day_Time'])
food_log=food_log.rename({'Food Name':'meal','id':'subject'},axis=1)
food_log['subject']=pd.to_numeric(food_log['subject'].str[-3:])
# use study meal or breakfast
studymeallist=['Study Meal','Study meal','study meal','Study']
nostudymeal_str=[]
for perid in food_log['subject'].unique():
    loctab=food_log[food_log['subject'].isin([perid])]
    flag_study=loctab['Group'].isin(studymeallist).any()
    if not flag_study:
        nostudymeal_str.append(perid)
study_meals=food_log[(food_log['Group'].isin(studymeallist)) & (~food_log['subject'].isin(nostudymeal_str))]
study_meals=pd.concat([study_meals, food_log[(food_log['Group'].isin(studymeallist+['Breakfast'])) & (food_log['subject'].isin(nostudymeal_str))]])
# np.sort(study_meals['subject'].unique())
# np.sort(food_log['subject'].unique())
# Clean cycles
study_meals['cycle'].unique()
study_meals.loc[study_meals['cycle'].isin(['cycle 1','Cycle 1']),'cycle']=1
study_meals.loc[study_meals['cycle'].isin(['cycle 2','Cycle 2', 'cycle 2.2','STUDYID-084 diet cycles']),'cycle']=2
study_meals.loc[study_meals['cycle'].isin(['cycle 3','Cycle 3']),'cycle']=3
study_meals=study_meals[study_meals['cycle']!=3]
# standarized meal names
rice = ["Trader Joe's, rice, jasmine, cooked, from frozen",
       'White Rice, Steamed', 'Jasmine Rice, Cooked in Unsalted Water',
       "Trader Joe's, Rice, Jasmine, Dry",
       "Trader Joe's, Organic Jasmine Rice",
       'White Rice, Cooked in Salted Water', 'study meal 3: rice',
       "Trader Joe's, rice, jasmine, Thai, cooked, from frozen"]
potatoes = ['Potatoes, frozen, whole, cooked, boiled, drained, without salt',
            "Trader Joe's, Shredded Hash Browns, Frozen",
           "Trader Joe's, Shredded Hash Browns",
           "Trader Joe's, Shredded Potato Hash Browns ",
           "Trader Joe's Shredded Potato Hash Browns",
           "Trader Joe's, Shredded Potato Hash Browns",
            'Potato sticks','Potato, Boiled without Skin','Potato, unknown preparation, without skin']
berries = ['Dole, strawberries, fresh', 'Veg-Fresh, Mixed Berries',
       'Blueberries, Fresh', 'Strawberries, Fresh',
       'Raspberries, Fresh, Red', 'Dole, Mixed Berries',
       'Blackberries, Fresh', "Driscoll's, Blackberries",
       'study meal 1: berries']
# bread = ['Oroweat, bread, country white', 'Oroweat, Country White Bread',
#        'White bread, store bought', 'Whole Wheat Bread, Store Bought',
#        'White Bread, Store Bought', 'Oroweat, Country Buttermilk Bread',
#        'study meal 4: bread', "Alfaro's, Artesano Style Bread",
#         'Orowheat, Thin-Sliced Rustic White']
bread = ['Oroweat, bread, country white', 'Oroweat, Country White Bread',
       'White bread, store bought', 
       'White Bread, Store Bought', 'Oroweat, Country Buttermilk Bread',
       'study meal 4: bread', "Alfaro's, Artesano Style Bread",
        'Orowheat, Thin-Sliced Rustic White']
pasta = ['Pasta, Dry, Enriched',
       'Angel Hair Pasta, Whole Wheat, Cooked in Unsalted Water',
        'Macaroni Noodles, White, Cooked in Unsalted Water',
        "Trader Joe's, Italian Macaroni","Pasta with extra fiber, cooked in unsalted water"]
grapes = ['Grapes, Fresh', 'study meal 2: grapes']
beans = ["Trader Joe's, Cuban Style Black Beans",
       'Black Beans, Canned, Drained',
       "Trader Joe's, black beans, cuban style, canned",
       "Trader Joe's, black beans, canned"]
fiber = ['Now, fiber, apple, powder', 'Equate, Fiber Powder', 'Pea fiber',
       'Leader, Fiber Powder', 'Fiber', 'Green Peas Fiber, Raw',
       "Barlean's, Fiber Blend, Acacia + Coconut + Flax",
       "Trader Joe's, Organic Pea Fiber Powder, Unsweetened",
       'Kroger, Fiber Powder ', 'Kroger, Fiber Powder',
       'Metamucil, Psyllium Fiber Supplement', 'Pea Fiber', 'study fiber',
       'Now Foods, Pea Fiber', 'Peak Bio, Max Fiber 2.0']
egg = ['Egg Whites Only, Cooked', 'egg whites', 'Boiled egg white',
       'study egg whites', 'Eggs, Cooked']
cream = ['Vermont Creamery, creme fraiche, madagascar vanilla',
       'Creme Fraiche','Vermont Creamery, creme fraiche, madagascar vanilla',
       'Whipping Cream, Extra Heavy, Gourmet, Not Whipped', 'cream',
       'Cream Cheese Spread',
       'Sour Cream', 'Isigny Ste Mre, Creame Fraiche', 'Study Cream',
       'study: cream', 'Whipped Cream, Sweetened',
       'Whipped Cream, Unsweetened']
glucose = ['Glucose Shot','Glucola']
quinoa = ['Quinoa, Cooked']
meal_names = rice+potatoes+berries+bread+pasta+grapes+beans+fiber+egg+cream+glucose+quinoa
standard_meals = {}
for i in rice:
    standard_meals[i] = 'rice'
for i in potatoes:
    standard_meals[i] = 'potatoes'
for i in berries:
    standard_meals[i] = 'berries'
for i in bread:
    standard_meals[i] = 'bread'
for i in pasta:
    standard_meals[i] = 'pasta'
for i in grapes:
    standard_meals[i] = 'grapes'
for i in beans:
    standard_meals[i] = 'beans'
for i in fiber:
    standard_meals[i] = 'fiber'
for i in egg:
    standard_meals[i] = 'egg'
for i in cream:
    standard_meals[i] = 'cream'
for i in glucose:
    standard_meals[i] = 'glucose'
for i in quinoa:
    standard_meals[i] = 'quinoa'
study_meals.loc[:,'meal']=study_meals['meal'].map(standard_meals)
# extract the main meal, not including mitigators
study_meals.loc[:,'Group']='study'
study_meals=study_meals.dropna(subset='meal')
mitigators=['egg','fiber','cream']
meal_no_mitigators=study_meals[~study_meals['meal'].isin(mitigators)]
meal_no_mitigators=meal_no_mitigators.sort_values(['subject','Day_Time'])
# meal_no_mitigators.groupby([pd.Grouper(key='Day_Time',freq='D'),'subject']).filter(lambda x: len(x)>1)
# All meals with more than 1 meal in a day (excluding mitigators) look good except for cycle 3. Which we will just ignore from now on.
# Reduce more than 1 meal down to just 1 (usually berries)
meal_no_mitigators=meal_no_mitigators.groupby([pd.Grouper(key='Day_Time', freq='D'),'subject']).head(1)
meal_no_mitigators['date']=meal_no_mitigators['Day_Time'].dt.date
study_meal_times=meal_no_mitigators[['Day_Time','subject', 'date']]
food_log['date']=food_log['Day_Time'].dt.date
# Food before the study meal, Labeling all as suspect but could probably keep some if they just have water
combined=pd.merge(study_meal_times,food_log,on=['subject','date'])
sus=combined[(~combined['meal'].isin(meal_names)) & (combined['Day_Time_x']>combined['Day_Time_y']) & (combined['Day_Time_y'].dt.hour>2)]
sus=sus[['subject','date']].value_counts().reset_index(name='count').drop('count',axis=1)
sus['food_before']=True
final_meals=pd.merge(meal_no_mitigators,sus,on=['subject','date'],how='left')
final_meals.loc[:,'food_before']=final_meals['food_before'].fillna(False)
# food within the range of studies CGM
combined=pd.merge(study_meal_times,food_log,on=['subject','date'])
sus=combined[(~combined['meal'].isin(meal_names)) & ((combined['Day_Time_y']-combined['Day_Time_x'])/np.timedelta64(1,'h')<=2) & (combined['Day_Time_y'].dt.hour>2) & (combined['Day_Time_y']>combined['Day_Time_x'])]
sus=sus[['subject','date']].value_counts().reset_index(name='count').drop('count',axis=1)
sus['food_within']=True
final_meals=pd.merge(final_meals,sus,on=['subject','date'],how='left')
final_meals.loc[:,'food_within']=final_meals['food_within'].fillna(False)
# Add mitigators
meal_mitigator=study_meals[study_meals['meal'].isin(mitigators)].copy()
meal_mitigator['date']=meal_mitigator['Day_Time'].dt.date
final_with_mitigators=pd.merge(final_meals,meal_mitigator[['Day_Time','meal','subject','date']],how='left',on=['subject','date'])
# remove incorrect mitigator meals (before carb, too much after carb)
flag_smallinterval=(final_with_mitigators['Day_Time_x']-final_with_mitigators['Day_Time_y'])/np.timedelta64(1,'h')<=0.5
flag_orderfood=(final_with_mitigators['Day_Time_x']-final_with_mitigators['Day_Time_y'])/np.timedelta64(1,'h')>=0
flag_nomitigator=final_with_mitigators['Day_Time_y'].isna()
final_with_mitigators=final_with_mitigators[(flag_smallinterval & flag_orderfood) | flag_nomitigator]
final_with_mitigators=final_with_mitigators.drop(['Day_Time_y'],axis=1).rename({'Day_Time_x':'Day_Time','meal_x':'meal','meal_y':'mitigator'},axis=1)
# 
final_with_mitigators['food']=final_with_mitigators['meal']+'+'+final_with_mitigators['mitigator'].fillna('')
final_with_mitigators['rep']=final_with_mitigators.groupby(['subject','food']).cumcount()+1
final_with_mitigators.to_pickle("meals_temp.pkl")
# collect information from food challenges
# food_log_cycle1=pd.read_csv(rawdatapath+'cronometer/meals.std.datetime_AfterPandemic.tsv',delimiter='\t',header=None,names=['id','time','food','rep'])
# food_log_cycle1['time']=pd.to_datetime(food_log_cycle1['time'])
# food_log_cycle1['subject']=pd.to_numeric(food_log_cycle1['id'].str[-3:])
# food_log_cycle1['rep']=pd.to_numeric(food_log_cycle1['rep'].str[-1:])
# food_log_cycle1=food_log_cycle1.drop('id',axis=1).rename({'time':'datetime_local'},axis=1)[['subject','datetime_local','food','rep']]
# food_log_cycle1['cycle']=1
# food_log_cycle1['foods']=food_log_cycle1['food']
# food_log_cycle1['food_before']=False
# food_log_cycle1['food_within']=False

new_meals=final_with_mitigators
new_meals=new_meals.rename({'Day_Time':'datetime_local'})
new_meals['meal']=new_meals['meal'].str.capitalize()
# new_meals=new_meals[(new_meals['cycle']==2) | (new_meals['cycle']==1)]
new_meals['mitigator']=new_meals['mitigator'].replace({'egg':'Protein','cream':'Fat','fiber':'Fiber'})
new_meals['food']=new_meals['meal']
new_meals['foods']=new_meals['meal']+('+'+new_meals['mitigator']).fillna('')
new_meals['datetime_local']=new_meals['Day_Time']
new_meals[['subject','datetime_local','food','foods','mitigator','rep','cycle']]
study_foods=pd.concat([new_meals[['subject','datetime_local','food','foods','mitigator','rep','cycle','food_before','food_within']]],axis=0,ignore_index=True)
# Add a couple that were missed in new processing
manual=pd.read_csv(proc_para_path+'manual_add_food.csv')
manual['datetime_local']=pd.to_datetime(manual['datetime_local'])
study_foods=pd.concat([study_foods,manual],ignore_index=True)
# modify some food record
timemodify=pd.read_csv(proc_para_path+'modify_time.tsv',sep='\t')
for index,row in timemodify.iterrows():
    ind_record=(study_foods['datetime_local'].dt.date==pd.to_datetime(row['dates']).date()) &(study_foods['subject']==row['subject'])
    study_foods.loc[ind_record,'datetime_local']=study_foods.loc[ind_record,'datetime_local']+pd.Timedelta(hours=row['shift'])
# 
# some subject specific change of time that's masked 
# 
study_foods=study_foods.sort_values(['subject','datetime_local'])
study_foods=study_foods.reset_index(drop=True)
study_foods['rep']=study_foods.groupby(['subject','foods']).cumcount()+1
study_foods.to_csv('study_foods.csv',index=False)
# merge with CGM data
cgm_df=pd.read_csv('final_cgm.csv',parse_dates=['datetime_local'])
study_foods=pd.read_csv('study_foods.csv',parse_dates=['datetime_local'])
new_dfs=[]
for index,row in study_foods.iterrows():
    start_meal_time=row['datetime_local']
    subject_df=cgm_df[cgm_df['subject']==row['subject']].dropna(subset='datetime_local').set_index('datetime_local')
    meal_df=subject_df.loc[start_meal_time-pd.Timedelta(minutes=90):start_meal_time+pd.Timedelta(minutes=180)].copy()
    meal_df['foods']=row['foods']
    meal_df['mitigator']=row['mitigator']
    meal_df['food']=row['food']
    meal_df['rep']=row['rep']
    meal_df = meal_df.reset_index()
    meal_df['mins_since_start']=(meal_df['datetime_local']-start_meal_time).dt.total_seconds()/60
    new_dfs.append(meal_df)
cgm_foods=pd.concat(new_dfs)
# Adjustments for daylight savings
# Find daylight savings times.  Look at any meals that occur within 10 days of a time switch.  If meals look off by 1 hour then switching the timing of the meals before joining to CGM and then switch it back after joined to CGM
# Identify daylight savings times
daylight=['2022-03-13','2022-11-06','2021-03-14','2021-11-07','2020-03-08','2020-11-01','2019-03-10','2019-11-03','2018-03-11','2018-11-04','2017-03-12','2017-11-05']
cgm_foods=cgm_foods.sort_values('datetime_local')
daylight_glucose=[]
for i in daylight:
    daylight_glucose.append(cgm_foods.set_index('datetime_local').loc[pd.to_datetime(i):pd.to_datetime(i)+pd.Timedelta(days=10)])
daylight_glucose=pd.concat(daylight_glucose)
daylight_glucose['subject'].unique()
daylight_meals=daylight_glucose.groupby(['subject','foods','rep']).size().reset_index()
# Load and adjust daylight savings
forward=pd.read_csv(proc_para_path+'/daylight_forward_modify.csv')
back=pd.read_csv(proc_para_path+'daylight_back.csv')
for index,row in forward.iterrows():
    food_challenge=study_foods[(study_foods['foods']==row['foods'])&(study_foods['rep']==row['rep'])&(study_foods['subject']==row['subject'])]
    if food_challenge.size>0:
        study_foods.at[food_challenge.index[0],'datetime_local']=study_foods.at[food_challenge.index[0],'datetime_local']+pd.Timedelta(hours=1)
for index,row in back.iterrows():
    food_challenge=study_foods[(study_foods['foods']==row['foods'])&(study_foods['rep']==row['rep'])&(study_foods['subject']==row['subject'])]
    if food_challenge.size>0:
        study_foods.at[food_challenge.index[0],'datetime_local']=study_foods.at[food_challenge.index[0],'datetime_local']-pd.Timedelta(hours=1)
# Join study foods and cgm together and revert time shifts to reflect local time
new_dfs=[]
for index,row in study_foods.iterrows():
    start_meal_time=row['datetime_local']
    subject_df=cgm_df[cgm_df['subject']==row['subject']].dropna(subset='datetime_local').set_index('datetime_local')
    meal_df=subject_df.loc[start_meal_time-pd.Timedelta(minutes=50):start_meal_time+pd.Timedelta(minutes=200)].copy()
    meal_df['foods']=row['foods']
    meal_df['mitigator']=row['mitigator']
    meal_df['food']=row['food']
    meal_df['rep']=row['rep']
    meal_df=meal_df.reset_index()
    meal_df['mins_since_start']=(meal_df['datetime_local']-start_meal_time).dt.total_seconds()/60
    new_dfs.append(meal_df)
cgm_foods=pd.concat(new_dfs)
# Looks good. Revert back start time of food challenges in case we want to use that information downstream
cgm_foods=cgm_foods.reset_index(drop=True)
for index,row in forward.iterrows():
    food_challenge=cgm_foods[(cgm_foods['foods']==row['foods'])&(cgm_foods['rep']==row['rep'])&(cgm_foods['subject']==row['subject'])]
    cgm_foods.loc[food_challenge.index,'datetime_local']=cgm_foods.loc[food_challenge.index,'datetime_local']-pd.Timedelta(hours=1)
for index,row in back.iterrows():
    food_challenge=cgm_foods[(cgm_foods['foods']==row['foods'])&(cgm_foods['rep']==row['rep'])&(cgm_foods['subject']==row['subject'])]
    cgm_foods.loc[food_challenge.index,'datetime_local']=cgm_foods.loc[food_challenge.index,'datetime_local']+pd.Timedelta(hours=1)
cgm_foods.to_csv('cgm_foods.csv',index=False)