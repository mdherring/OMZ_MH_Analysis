# 00_functions

# This script includes all of the functions for this project.

# create_df_cols(): if one of the dramv databases (kegg, viral, pfam, vogdb) doesn't have a column, create one with 0 as the value
def create_db_cols(df): ''' input = a dataframe '''
    for col in df:
        if 'kegg_hit' not in df.columns:
            df.insert(loc=len(df.columns), column='kegg_hit', value=0)
        if 'viral_hit' not in df.columns:
            df.insert(loc=len(source_df.columns), column='viral_hit', value=0)
        if 'pfam_hits' not in df.columns:
            df.insert(loc=len(source_df.columns), column='pfam_hits', value=0)
        if 'vogdb_hits' not in df.columns:
            df.insert(loc=len(source_df.columns), column='vogdb_hits', value=0)
    return df

# split_classification(): split a column named 'classification' into taxonomy designations and adds them to the data data columns formatted without an underscore
def split_classification(df): ''' input = a data frame '''
    df[['domain', 'phyla', 'class', 'order', 'family', 'genus', 'species']] = df.classification.str.split(';', expand=True)
    df['domain'] = df['domain'].str.replace('d__', '') # remove the d__ in front of all observations
    df['phyla'] = df['phyla'].str.replace('p__', '')
    df['class'] = df['class'].str.replace('c__', '')
    df['order'] = df['order'].str.replace('o__', '')
    df['family'] = df['family'].str.replace('f__', '')
    df['genus'] = df['genus'].str.replace('g__', '')
    df['species'] = df['species'].str.replace('s__', '')
    return df

# create_counts_row(): create a dataframe with each row a different fasta and the columns for the number of genes annotated by each database as well as a total number of genes annotated column
def create_count_row(file): ''' input = a file '''
    df = pd.read_csv(file)
    df_name = df.iloc[0, 0]
    df_count = pd.DataFrame(df['annotation_source'].value_counts().reset_index())
    df_count['ID'] = df_name
    df_piv = df_count.pivot(index='ID', columns='index', values='annotation_source')
    df_piv = create_db_cols(df_piv)
    df_piv.rename(columns={'kegg_hit':'kegg_count','pfam_hits': 'pfam_count','vogdb_hits':'vogdb_count'}, inplace=True)
    df_piv['total_genes_annot'] = df_piv['kegg_count'] + df_piv['viral_hit'] + df_piv['pfam_count'] + df_piv['vogdb_count']
    df_piv['V_count'] = len(df[df['V'] == 1])
    df_piv['M_count'] = len(df[df['M'] == 1])
    df_piv['A_count'] = len(df[df['A'] == 1])
    df_piv['P_count'] = len(df[df['P'] == 1])
    df_piv['E_count'] = len(df[df['E'] == 1])
    df_piv['K_count'] = len(df[df['K'] == 1])
    df_piv['T_count'] = len(df[df['T'] == 1])
    df_piv['F_count'] = len(df[df['F'] == 1])
    df_piv['B_count'] = len(df[df['B'] == 1])
    return df_piv

# count_concat(): apply create_count_row() to multiple files using a pattern and concatenate the results into one data frame
def count_concat(csv_pattern): ''' input = file pattern '''
    csv_file_paths = glob.glob(csv_pattern)
    dfs_list = []
    for file in csv_file_paths:
        df_piv = create_count_row(file)
        dfs_list.append(df_piv)
    result_df = pd.concat(dfs_list)
    return result_df