# 13_clean_dramv_annot_script

# This script applies the functions for the jupyter notebook dramv_function_tests. This code prepares all dramv annotation output files for summarizing and formats them uniformly.

import pandas as pd
import math
import glob
import os # these two packages are good for searching and navigating file systems
import os.path as op

# get_ann_text function:
def get_ann_text(hit_text, column_type = 'viral_hit'):
    if type(hit_text) == float:
        return hit_text
    if column_type == 'viral_hit':
        no_org = hit_text.split("[")[0]
        no_acc_id = " ".join(no_org.split(" ")[1:-1]) 
        return no_acc_id
    if column_type in ['kegg_hit']:
        no_ee = hit_text.split("[")[0].strip()
        return no_ee
    if column_type == 'pfam_hits':
        no_pf_ids = ";".join([text.split("[")[0].strip() for text in hit_text.split(";")]) 
        return no_pf_ids
    if column_type == 'vogdb_hits':
        no_code = hit_text.split(";")[0]
        no_acc = " ".join(no_code.split(" ")[1:])
        return no_acc
    
# assign_annot function:
def assign_annot(line):
        annot_source = math.nan
        if line['rank'] == 'A' :
            annot_source = 'kegg_hit'
        elif line['rank'] == 'B' :
            annot_source = 'viral_hit'
        elif line['rank'] == 'C' and type(line['kegg_hit']) == str :
            annot_source = 'kegg_hit'
        elif line['rank'] == 'C' and not pd.isna(line['kegg_hit']) and type(line['viral_hit']) == str :
            annot_source = 'viral_hit'
        elif line['rank'] == 'C' and pd.isna(line['kegg_hit']) and type(line['viral_hit']) == str :
            annot_source = 'vogdb_hits'
        elif line['rank'] == 'D' :
            annot_source = 'pfam_hits'
        elif line['rank'] == 'E' and not pd.isna(line['vogdb_hits']) and type(line['vogdb_hits']) == str :
            annot_source = 'vogdb_hits'                          
        elif line['rank'] == 'E' and pd.isna(line['vogdb_hits']) and type(line['vogdb_hits']) == str :
            annot_source = np.na
        else:
            return math.nan, math.nan
        keep_annot = get_ann_text(line[annot_source], column_type = annot_source)
        return keep_annot, annot_source
    
# convert_na_to_float function:
def convert_na_to_float(df):
    for col in df.columns:
        for i in range(len(df)):
            if pd.isna(df.at[i, col]):
                df.at[i, col] = float('NaN')

# for loop: create database columns if they don't already exist and put NAs then convert NAs to float type for get_ann_text function, apply get_ann_text and assign_annot functions, create a file with only certain columns
tsv_pattern = "/Users/melissaherring/Google Drive/My Drive/MH_project/dramv/*/annotations.tsv"  # Replace with your file pattern
tsv_file_paths = glob.glob(tsv_pattern)
columns_to_keep = ['fasta', 'scaffold', 'start_position', 'end_position', 'annotation','annotation_source','amg_flags','V', 'M', 'A', 'P', 'E', 'K', 'T', 'F', 'B']
num_columns_list = []


for file_path in tsv_file_paths:
    df = pd.read_csv(file_path, delimiter='\t')
    
    if 'kegg_hit' not in df.columns:
        df.insert(loc=len(df.columns), column='kegg_hit', value=pd.NA)
    if 'viral_hit' not in df.columns:
        df.insert(loc=len(df.columns), column='viral_hit', value=pd.NA)
    if 'pfam_hits' not in df.columns:
        df.insert(loc=len(df.columns), column='pfam_hits', value=pd.NA)
    if 'vogdb_hits' not in df.columns:
        df.insert(loc=len(df.columns), column='vogdb_hits', value=pd.NA)
        
    convert_na_to_float(df)
            
    df['viral_ann_text'] = df['viral_hit'].apply(get_ann_text, args = ('viral_hit',))
    df['kegg_ann_text'] = df['kegg_hit'].apply(get_ann_text, args = ('kegg_hit',)) 
    df['pfam_ann_text'] = df['pfam_hits'].apply(get_ann_text, args = ('pfam_hits',)) 
    df['vogdb_ann_text'] = df['vogdb_hits'].apply(get_ann_text, args = ('vogdb_hits',))
    
    df[['annotation','annotation_source']] = df.apply(assign_annot, axis=1, result_type='expand')
    
    df = df.applymap(str)

    new_column_names = list(['V', 'M', 'A', 'P', 'E', 'K', 'T', 'F', 'B'])

    flag_columns = defaultdict(lambda: [])

    for i in new_column_names:
        for flag in df['amg_flags']:
            
            if i in flag:
                flag_columns[i].append(1)
            else:
                flag_columns[i].append(0)

    newdf = pd.concat([df, pd.DataFrame(flag_columns)], axis = 1)
    
    newdf_name = newdf.iloc[0, 1]
    dir_path = '/Users/melissaherring/Google Drive/My Drive/MH_project/dramv_trim/'
    file_name = f"{newdf_name}.csv"
    full_path = dir_path + file_name
    newdf[columns_to_keep].to_csv(full_path, index = False)