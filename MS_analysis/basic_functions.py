import pandas as pd
import numpy as np


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
