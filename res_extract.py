import pandas as pd
import numpy as np

pd.set_option('display.max_rows', None)
pd.options.mode.chained_assignment = None
drop_cols = [5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22]
keep_cols = [
    '#CHROM',
    'POS',
    'REF',
    'ALT',
    '#CHROM_2',
    '#CHROM_7',
]

def extract_res(clean_up):
    split_df = pd.concat(
        [clean_up['INFO'].str.split(
            ';', expand = True
        ) for x in ['INFO']], axis = 1, keys = clean_up.columns
    )
    split_df.columns = split_df.columns.map(
        lambda x: '_'.join((x[0], str(x[1])))
    )
    split_df = split_df.replace({'':np.nan, None:np.nan})
    cat_df = pd.concat([clean_up, split_df], axis = 1)
    reindex_cat_df = cat_df[keep_cols] # can bring this lower for better clean up
    reindex_cat_df['REFPOSALT'] = reindex_cat_df['REF'] + reindex_cat_df['POS'].astype(str) + reindex_cat_df['ALT']
    split_again = pd.concat(
        [reindex_cat_df['#CHROM_7'].str.split('|', expand = True
    ) for x in ['#CHROM_7']], axis = 1, keys = reindex_cat_df.columns)
    split_again.columns = split_again.columns.map(
        lambda x: '_'.join((x[0], str(x[1])))
    )
    re_split_again = split_again.replace({'':np.nan, None:np.nan})
    re_cat_df = pd.concat([reindex_cat_df, re_split_again], axis = 1)
    re_cat_df.drop(re_cat_df.columns[drop_cols], axis = 1, inplace = True)
    re_cat_df.rename(columns = {'#CHROM_3':'Locus_Tag', '#CHROM_10':'PROT-POS'}, inplace = True) #column name changes
    re_cat_df['PROT-POS'] = re_cat_df['PROT-POS'].map(lambda x: x.lstrip('p.'))
    return re_cat_df