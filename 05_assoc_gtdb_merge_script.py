# 05_assoc_gtdb_merge_script

# This python script merges the MAG GTDB data with the host-virus associations. The resulting data frame consists of each host-virus association (rows) with all of the columns of the host-virus associations along with their phylogenic classifications.

# Load packages
import pandas as pd
import os
import sys
import csv
import numpy as np
import matplotlib
import glob
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt

# Load data
associations = pd.read_csv('~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/all_host_associations.csv')
mags = jv154_bac = pd.read_csv('~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/mag_data/all_mag_gtdb.csv')

# rename MAGs user_genome column to match associations data frame
mags.rename(columns={'user_genome':'host_name'}, inplace=True) 

# Merge two data frames together
combo = associations.merge(mags,how='left',on=["host_name","sample_name","sample_depth"])

# make all NA values match
combo = combo.fillna("NA")

# write csv
combo.to_csv('~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/all_associations_gtdb.csv', index=False)