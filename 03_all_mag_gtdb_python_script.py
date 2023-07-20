# 03_all_mag_gtdb_python_script

# Run this python script to create a merged data frame of both the arachaea and bacteria MAGs. Data is original files from google drive.

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

jv119_arc = pd.read_csv('./Data/mag_data/gtdbtk/jv-119/jv-119.ar53.summary.tsv',sep='\t')
jv119_bac = pd.read_csv('./Data/mag_data/gtdbtk/jv-119/jv-119.bac120.summary.tsv',sep='\t')
jv121_arc = pd.read_csv('./Data/mag_data/gtdbtk/jv-121/jv-121.ar53.summary.tsv',sep='\t')
jv121_bac = pd.read_csv('./Data/mag_data/gtdbtk/jv-121/jv-121.bac120.summary.tsv',sep='\t')
jv132_arc = pd.read_csv('./Data/mag_data/gtdbtk/jv-132/jv-132.ar53.summary.tsv',sep='\t')
jv132_bac = pd.read_csv('./Data/mag_data/gtdbtk/jv-132/jv-132.bac120.summary.tsv',sep='\t')
jv154_arc = pd.read_csv('./Data/mag_data/gtdbtk/jv-154/jv-154.ar53.summary.tsv',sep='\t')
jv154_bac = pd.read_csv('./Data/mag_data/gtdbtk/jv-154/jv-154.bac120.summary.tsv',sep='\t')

# create sample_name column
jv119_arc['sample_name'] = "JV119"
jv119_bac['sample_name'] = "JV119"
jv121_arc['sample_name'] = "JV121"
jv121_bac['sample_name'] = "JV121"
jv132_arc['sample_name'] = "JV132"
jv132_bac['sample_name'] = "JV132"
jv154_arc['sample_name'] = "JV154"
jv154_bac['sample_name'] = "JV154"

# create sample_depth column
jv119_arc['sample_depth'] = 400
jv119_bac['sample_depth'] = 400
jv121_arc['sample_depth'] = 95
jv121_bac['sample_depth'] = 95
jv132_arc['sample_depth'] = 80
jv132_bac['sample_depth'] = 80
jv154_arc['sample_depth'] = 140
jv154_bac['sample_depth'] = 140

# combine all data frames into one data set
combo = pd.concat([jv119_arc, jv119_bac, jv121_arc, jv121_bac, jv132_arc, jv132_bac, jv154_arc, jv154_bac])

# split the classification into tax level columns and rename
combo[['domain', 'phyla', 'class', 'order', 'family', 'genus', 'species']]=combo.classification.str.split(';', expand=True)
combo['domain'] = combo['domain'].str.replace('d__', '') # remove the d__ in front of all observations
combo['phyla'] = combo['phyla'].str.replace('p__', '')
combo['class'] = combo['class'].str.replace('c__', '')
combo['order'] = combo['order'].str.replace('o__', '')
combo['family'] = combo['family'].str.replace('f__', '')
combo['genus'] = combo['genus'].str.replace('g__', '')
combo['species'] = combo['species'].str.replace('s__', '')