# 08_vMAG_associations_script

import pandas as pd
import os
import sys
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import seaborn as sns
from collections import Counter

assoc = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/all_associations_gtdb.csv')
assoc_vMAG = assoc[assoc['virus_type'] == 'vMAG']

vMAG_119 = pd.read_csv("~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-119_874814_Viral_Files/viral_MAGs/viral_mags_summary.tsv", sep = '\t')
vMAG_121 = pd.read_csv("~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-121_874818_Viral_Files/viral_MAGs/viral_mags_summary.tsv", sep = '\t')
vMAG_132 = pd.read_csv("~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-132_874826_Viral_Files/viral_MAGs/viral_mags_summary.tsv", sep = '\t')
vMAG_154 = pd.read_csv("~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-154_874822_Viral_Files/viral_MAGs/viral_mags_summary.tsv", sep = '\t')

# combine all vMAG data frames together
vMAGs_all = pd.concat([vMAG_119,vMAG_121,vMAG_132,vMAG_154])

# split up contig_id column into different columns
vMAGs_all['virus_name'] = vMAGs_all['contig_id'].str.split("|").str[0]
vMAGs_all['N'] = vMAGs_all['contig_id'].str.split("|").str[1]
vMAGs_all['N'] = vMAGs_all['N'].str.replace('N=', '')
vMAGs_all['N'] = pd.to_numeric(vMAGs_all['N'])
vMAGs_all['L'] = vMAGs_all['contig_id'].str.split("|").str[2]
vMAGs_all['L'] = vMAGs_all['L'].str.replace('L=', '')
vMAGs_all['L'] = pd.to_numeric(vMAGs_all['L'])

# merge vMAGs with associations 
combo = assoc_vMAG.merge(vMAGs_all,how='inner',on=['virus_name'])

combo.to_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/vMAG_associations.csv', index=False)