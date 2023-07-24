# all_host_virus_association_python_script

# This script creates a csv file that consists of all 4 sample's contigs and vMAGs combined into 1 dataset.

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

jv119_contigs = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-119_874814_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv', sep ='\t')
jv121_contigs = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-121_874818_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv', sep ='\t')
jv132_contigs = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-132_874826_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv', sep ='\t')
jv154_contigs = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-154_874822_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv', sep = '\t')

jv119_vMAG = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-119_874814_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv', sep ='\t')
jv121_vMAG = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-121_874818_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv', sep ='\t')
jv132_vMAG = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-132_874826_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv', sep ='\t')
jv154_vMAG = pd.read_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-154_874822_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv', sep ='\t')

# Create sample_name column for each data frame

jv119_contigs['sample_name'] = 'JV119'
jv121_contigs['sample_name'] = 'JV121'
jv132_contigs['sample_name'] = 'JV132'
jv154_contigs['sample_name'] = 'JV154'

jv119_vMAG['sample_name'] = 'JV119'
jv121_vMAG['sample_name'] = 'JV121'
jv132_vMAG['sample_name'] = 'JV132'
jv154_vMAG['sample_name'] = 'JV154'

# Create virus_type column for each data frame

jv119_contigs['virus_type'] = 'contig'
jv121_contigs['virus_type'] = 'contig'
jv132_contigs['virus_type'] = 'contig'
jv154_contigs['virus_type'] = 'contig'

jv119_vMAG['virus_type'] = 'vMAG'
jv121_vMAG['virus_type'] = 'vMAG'
jv132_vMAG['virus_type'] = 'vMAG'
jv154_vMAG['virus_type'] = 'vMAG'

# Create sample_depth column for each data frame

jv119_contigs['sample_depth'] = 400
jv119_vMAG['sample_depth'] = 400
jv121_contigs['sample_depth'] = 95
jv121_vMAG['sample_depth'] = 95
jv132_contigs['sample_depth'] = 80
jv132_vMAG['sample_depth'] = 80
jv154_contigs['sample_depth'] = 140
jv154_vMAG['sample_depth'] = 140

# rename columns to match

df_list1 = [jv119_contigs, jv121_contigs, jv132_contigs, jv154_contigs]

column_mapping_1 = {
    'mobile_contig_name': 'virus_name',
    'mobile_contig_length (bp)': 'virus_length',
    'mobile_contig_read_count (reads)': 'virus_read_count',
    'mobile_contig_read_depth (reads/kbp)': 'virus_read_depth',
    'mobile_contig_read_depth_in_this_cluster (reads/kbp)': 'virus_read_depth_in_host',
    'cluster_name': 'host_name',
    'cluster_length (bp)': 'host_length',
    'cluster_read_count (reads)': 'host_read_count',
    'cluster_read_depth (reads/kbp)': 'host_read_depth',
    'intra_read_count (reads)': 'intra_read_count',
    'intra_linkage_density (reads/kbp^2)': 'intra_linkage_density',
    'inter_read_count (reads)': 'inter_read_count',
    'raw_inter_linkage_density (reads/kbp^2)': 'raw_inter_linkage_density',
    'raw_inter_vs_intra_ratio': 'raw_inter_vs_intra_ratio',
    'mobile_element_copies_per_cell': 'viral_copies_per_cell',
    'adjusted_inter_connective_linkage_density (reads/kbp^2)': 'adjusted_inter_linkage_density',
    'adjusted_inter_vs_intra_ratio': 'adjusted_inter_vs_intra_ratio',
    'sample_name':'sample_name',
    'virus_type': 'virus type',
    'sample_depth': 'sample_depth'
}

for df in df_list1:
    df.rename(columns=column_mapping_1, inplace=True)
    
df_list2 = [jv119_vMAG, jv121_vMAG, jv132_vMAG, jv154_vMAG]

column_mapping_2 = {
    'mobile_cluster_name': 'virus_name',
    'mobile_cluster_length (bp)': 'virus_length',
    'mobile_cluster_read_count (reads)': 'virus_read_count',
    'mobile_cluster_read_depth (reads/kbp)': 'virus_read_depth',
    'mobile_cluster_read_depth_in_this_cluster (reads/kbp)': 'virus_read_depth_in_host',
    'cluster_name': 'host_name',
    'cluster_length (bp)': 'host_length',
    'cluster_read_count (reads)': 'host_read_count',
    'cluster_read_depth (reads/kbp)': 'host_read_depth',
    'intra_read_count (reads)': 'intra_read_count',
    'intra_linkage_density (reads/kbp^2)': 'intra_linkage_density',
    'inter_read_count (reads)': 'inter_read_count',
    'raw_inter_linkage_density (reads/kbp^2)': 'raw_inter_linkage_density',
    'raw_inter_vs_intra_ratio': 'raw_inter_vs_intra_ratio',
    'mobile_element_copies_per_cell': 'viral_copies_per_cell',
    'adjusted_inter_connective_linkage_density (reads/kbp^2)': 'adjusted_inter_linkage_density',
    'adjusted_inter_vs_intra_ratio': 'adjusted_inter_vs_intra_ratio',
    'sample_name':'sample_name',
    'virus_type': 'virus type',
    'sample_depth': 'sample_depth'
}

for df in df_list2:
    df.rename(columns=column_mapping_2, inplace=True)

# merge data frames
combo = pd.concat([jv119_contigs, jv121_contigs, jv132_contigs, jv154_contigs, jv119_vMAG, jv121_vMAG, jv132_vMAG, jv154_vMAG])

# write csv
#combo.to_csv('~/Documents/Bigelow/Virus_Project/OMZ_MH_Analysis/Data/proximeta_viral_files/all_host_associations.csv', index=False)