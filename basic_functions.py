import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import pairwise2

def read_input_filter_proteotypic(df):
    '''''
    reads a csv file with MS data (for example format see LIP school Report). Filters out all non-proteotypic peptides. Drops duplicates based on condition, replicate and precursorID. Returns dataframe with the important columns.
    USAGE: df = read_input_filter_proteotypic("FILENAME OR VARIABLE")
    '''''
    df = pd.read_csv(df)
    df_filtered = df.loc[(df['PEP.IsProteotypic'] != False)]
    df_filtered = df_filtered.drop_duplicates(["R.Condition","R.Replicate","EG.PrecursorId"])
    df_filtered = df_filtered[["R.Condition","R.Replicate","PG.ProteinAccessions","PEP.StrippedSequence","FG.Quantity","EG.PrecursorId"]]
    return df_filtered




def split_conditions(df,cond1,cond2):
    '''''
    splits dataframe according to experimental condition. returns two dataframes, one per condition.
    USAGE example: df1, df2 = split_conditions(df2,"control","rapamycin")
    '''''
    df_cond1 = df[df["R.Condition"]==cond1]
    df_cond2 = df[df["R.Condition"]==cond2]
    return df_cond1,df_cond2



def average_replicates(dfin,merge=False,keep=True):
    '''''
    groups dataframe by precursorID, then computes the mean and stddev of the column FG.Quantity. if merge is set to TRUE, the mean and stddev are added to the input dataframe.
    USAGE: df = average_replicates(dfin,merge=True/False) (default = False)
    '''''
    df = dfin["FG.Quantity"].groupby(dfin["EG.PrecursorId"])
    df_ = df.mean()
    df__ = df.std()
    df_average = pd.DataFrame(df_)
    df_average["STDDEV"] = pd.DataFrame(df__)
    if merge:
        dfm = df_average.merge(dfin,on="EG.PrecursorId")
        if keep:
            return dfm
        else:
            dfm2 = dfm.drop_duplicates(["EG.PrecursorId","FG.Quantity_x","STDDEV","R.Condition"])
            return dfm2
    else:
        return df_average

def append_value(dict_obj, key, value):
    '''''
    This function is useful for creating a dicitonary where new values can be added and do not replace the value if the key was already in the dict.
    Use instead of dict.update()
    USAGE: append_value(dict,key,value). Can also be as dict inside a dict: append_value(pep_coord,k,{x:pep_coordinates})

    '''''
    # Check if key exist in dict or not
    if key in dict_obj:
        # Key exist in dict.
        # Check if type of value of key is list or not
        if not isinstance(dict_obj[key], list):
            # If type is not list then make it list
            dict_obj[key] = [dict_obj[key]]
        # Append the value in list
        dict_obj[key].append(value)
    else:
        # As key is not in dict,
        # so, add key-value pair
        dict_obj[key] = value

def align_pep_upseq(pep,upseq):
    alignments = pairwise2.align.localms(pep, upseq,1, -1, -1, 0)
    score = alignments[0][2]
    start = alignments[0][3]
    stop = alignments[0][4]
    return [score,start,stop]
