import datetime
from io import StringIO
import pandas as pd
import res_extract as re

def file_cleanup(file):
    """
    Remove the lines of the vcf file that contain the ##
    """
    with open(file, 'r') as vcf:
        oneline = ''
        lines = vcf.readlines()
        for line in lines:
            if not line.startswith('##'):
                oneline += line
        return oneline

def file2df(file):
    """
    Sets up the read of vcf file without the seriously unnecessary hashes,
    also will print which file is being read for the log.
    """
    now = datetime.datetime.now()
    time_log = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"{time_log}: Reading File - {file}")
    return pd.read_csv(StringIO(file_cleanup(file)), sep='\t', header = 0)

def reference_split(file):
    res_df = pd.DataFrame(file2df(file))
    if res_df.empty:
        return res_df
    if res_df["#CHROM"].iloc[0] == "FQ312006":
        return for_FQ(file)
    elif res_df["#CHROM"].iloc[0] == "":
        return for_L(file)
    else:
        raise Exception ("Incompatible Reference")

def for_FQ(file):
    res_df = pd.DataFrame(file2df(file))
    filt_res_df = res_df[res_df['INFO'].str.contains("HIB_12890|HIB_14190|HIB_16880|HIB_16890")]
    clean_up = filt_res_df[~filt_res_df['INFO'].str.contains("intergenic|synonymous")]
    if clean_up.empty:
        return clean_up
    return re.extract_res(clean_up)

def for_L(file):
    res_df = pd.DataFrame(file2df(file))
    filt_res_df = res_df[res_df['INFO'].str.contains("HI_1132|HI_1264|HI_1529|HI_1528")]
    clean_up = filt_res_df[~filt_res_df['INFO'].str.contains("intergenic|synonymous")]
    if clean_up.empty:
        return clean_up
    return re.extract_res(clean_up)
