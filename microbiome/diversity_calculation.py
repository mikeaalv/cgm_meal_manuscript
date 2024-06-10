# Adapated From Andrew Wallace Brooks' jupternotebook 
# calculate diversity alpha and beta for the microbiome data
from IPython.display import display, HTML # Jupyter Environment Display Options
from IPython.core.interactiveshell import InteractiveShell # Evaulate All
InteractiveShell.ast_node_interactivity = "all"
from IPython.utils import io # Capture the output of cells and write to file
import warnings # Used to Capture an Supress Warnings
from IPython.display import clear_output
# 
import os # Tool for terminal and operating system type calls
import sys # System tools for import variables
import glob # Tool to Regex Search for Files
import itertools # Iterate through data
import shutil # Remove directory
import time # Time and Date Tools
import random # Generate Random Values
import re
import string # Manipulate string objects
import copy # Generate copy and deepcopy of objects
from decimal import Decimal # Replacement for float when precision is necessary
from typing import OrderedDict # Dictionary that maintains key order
import math # Basic mathematics functions
from collections.abc import Iterable # Check if variable is iterable
import subprocess # Tool to make terminal calls from python
random.seed(19) # Set Random Seed for Reproducibility
# 
import numpy as np # Numpy Numerical Toolkit
import pandas as pd # Pandas Dataframes
import scipy as sp  # Scipy Scientific Toolkit
import patsy # Pandas Regression Formatting
import sklearn as sk
# 
import matplotlib.pyplot as plt # Main python plotting package interface
import seaborn as sns # Seaborn advanced dataframe plotting toolkit
import skbio
# Show figures within notebook instead of opening new window
%matplotlib inline
### CHANGE PANDAS DEFAULT DISPLAY FUNCTIONALITY
### COLUMNS ###
pd.set_option('display.max_columns',None)
pd.set_option('max_colwidth',400)
### ROWS ###
pd.set_option('display.max_rows',20)
### USE CSS TO SET NOTEBOOK DISPLAY STYLES ###
display(HTML(data="""
<style>
    /* ### SET NOTEBOOK TO FILL BROWSER WIDTH ### */
    div#notebook-container    { width: 100%; }
    div#menubar-container     { width: 100%; }
    div#maintoolbar-container { width: 100%; }
</style>
"""))
### PRINT DATE AND TIME ###
print('Date: '+time.strftime("%d/%m/%Y")+' '+time.strftime("%H:%M:%S"))
print('All Packages Loaded and Environment Initialized...')
### CURRENT COUNTS TABLE ###
pardir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/result/"
### OUTPUT FOLDER ###
outfolder=pardir+"/microbiome/diversity/"
# Import table #
tablepy=pd.read_csv(pardir+"microbiome/matrix_res/species_full.tsv",sep='\t',header=1)
### BIOM AND QIIME NEED:
# 1. table with columns as cleaner (shorter) sample names. And taxonomy as the last column
# 2. for each sample to add to 1 (i.e. each taxonomy is a proportion out of a total of 1 for each sample)
# 3. The third issue is each taxonomic level is the sum of itself and all more specific taxonomies. We only want each taxonomy to have the counts directly assigned to it.
# Rename Taxonomy Columns #
tablepy=tablepy.rename({'clade_name':'taxonomy','NCBI_tax_id':'tax_id'},axis=1)
# Re-Move Taxonomy #
tablecounts=tablepy.iloc[:,2:]# + tablepy[['taxonomy','tax_id']]
# Shorten Column Names
# tablecounts.columns=['-'.join(col2.split('-')[0:-1]) for col2 in ['-'.join(col.split('_')[:2]) if len(col.split('_'))==6 else col.split('_')[:1][0] for col in tablecounts.columns]]
# Get Sample Names #
samplecols=tablecounts.columns
# Re-Add Taxonomy #
tablecounts['taxonomy']=tablepy['taxonomy']
# Taxonomic Indexes #
taxindexes={}
for idx,tax in tablecounts['taxonomy'].items(): taxindexes[tax]=idx
### Now convert to a total of 100% for each sample by subtracting % of more specific assignments from less specific ones ###
sampledict={}
# For each Sample #
for cursam in samplecols:
  sampledict[cursam]={}
  # For each taxonomy (in reverse order, more specific to less specific assignments) #
  for curidx, curtax in tablecounts['taxonomy'].iloc[::-1].items():
    subtract=0
    # Loop through taxa already added and subtract abundance if they are a more specific assignment
    for prevtax in sampledict[cursam].keys():
      if curtax+"|" in prevtax:
        subtract -= sampledict[cursam][prevtax]
    sampledict[cursam][curtax]=tablecounts.loc[curidx,cursam]+subtract
    if sampledict[cursam][curtax] < 0:
       sampledict[cursam][curtax]=0
# Now all of the percents add up to 100% and the most accurate taxonomy is the only one with the counts #
finalotus=pd.DataFrame(sampledict)
# Lets convert to proportions (0-1):
finalotus=finalotus/finalotus.sum()
# Reset taxonomy index to column and name taxonomy #
finalotus=finalotus.reset_index().rename({'index':'taxonomy'},axis=1)
# Store Taxonomy #
curtax=finalotus['taxonomy']
# Remove taxonomy #
finalotus=finalotus.iloc[:,1:]
# Add taxonomy as last column for biom program #
finalotus['taxonomy']=curtax
# Save reformated matrix
finalotus.to_csv(outfolder+'otu_table.tsv',sep='\t',index_label='index',float_format='%.30f')
# Get counts as array
otuarray=finalotus.drop('taxonomy',axis=1).T.values
otucols=finalotus.drop('taxonomy',axis=1).columns
# Create Dictionaries to store different metrics #
alpha_diversity={}
beta_diversity={}
### Alpha Diversity ###
for alpha_metric in ['chao1','gini_index','observed_otus','shannon','simpson']:
    alpha_diversity[alpha_metric]=skbio.diversity.alpha_diversity(alpha_metric,otuarray,otucols)
### Beta Diversity ###
beta_diversity['braycurtis']=skbio.diversity.beta_diversity('braycurtis',otuarray,otucols)
beta_diversity['jaccard']=skbio.diversity.beta_diversity('jaccard',otuarray,otucols)
# Lets check how correlated our beta diversity is between weighted and unweighted metrics #
r,p_value,n=skbio.stats.distance.mantel(beta_diversity['braycurtis'],beta_diversity['jaccard'])
print(r)
print(p_value)
### PCOA - Like PCA ###
pcoa_bc=skbio.stats.ordination.pcoa(beta_diversity['braycurtis'])
print(pcoa_bc.proportion_explained)
pcoa_jaccard=skbio.stats.ordination.pcoa(beta_diversity['jaccard'])
print(pcoa_jaccard.proportion_explained)
# 
bc_pcoa_df=pcoa_bc.samples.iloc[:,0:5]
bc_pcoa_df=bc_pcoa_df.add_prefix('bc')
jac_pcoa_df=pcoa_jaccard.samples.iloc[:,0:5]
jac_pcoa_df=jac_pcoa_df.add_prefix('jaccard')
div_tab=pd.concat([pd.DataFrame(alpha_diversity),bc_pcoa_df,jac_pcoa_df],axis=1)
div_tab.to_csv(outfolder+'diversity_table.tsv',sep='\t',index_label='index',float_format='%.30f')