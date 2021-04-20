import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import pairwise2
import re
import json


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

def compute_sequence_coverage(df,accession="PG.ProteinAccessions_x",peptide_seq="PEP.StrippedSequence",merge=True):
    '''''
    computes sequence coverage from a dataframe. give the column name of the UP accession numbers and the column name of the peptide sequence.

    '''''
    with open("/Users/fsesterhenn/Dropbox/PostDoc/databases/UNIPROT_functional_annotations/UPID_sequences.json") as fp:
        UPseq  = json.load(fp)
    acc_pep = convert_df_to_dict(df,accession,peptide_seq)
    out = {}
    for k,v in acc_pep.items():
        try:
            total_length = len(UPseq[k])
            ar = np.zeros(total_length)
        except:
            continue
        for p in v:
            m = re.search(p,UPseq[k])
            ar[m.start():m.end()]=1
        coverage = sum(ar)/len(ar)
        out.update({k:coverage})
    if merge:
        df["sequence_coverage"] = df[accession].map(out)
        return df
    else:
        return out


def istryptic(df,column_name="PEP.StrippedSequence",ids = "PG.ProteinAccessions"):
    with open("/Users/fsesterhenn/Dropbox/PostDoc/databases/UNIPROT_functional_annotations/UPID_sequences.json") as fp:
        UPseq  = json.load(fp)
    istryptic = []
    peps = df[column_name].to_list()
    ids = df[ids].to_list()
    for i,p in enumerate(peps):
        length_begin = len(istryptic)
        try:
            fullseq = UPseq[ids[i]]

            match = re.search(p,fullseq)
            start = match.start()
            end = match.end()
            if p[-1] in ["K","R"]:
                if fullseq[start-1] in ["K","R"]:
                    istryptic.append('tryptic')
                else:
                    if start < 2:
                        istryptic.append("Prot_N_terminal")
                    else:
                        istryptic.append('PK_Cterm')
            elif p[-1] not in ["K","R"]:
                if end/len(fullseq) > 0.95:
                    istryptic.append("Prot_C_terminal")
                elif fullseq[start-1] in ["K","R"]:
                    istryptic.append('PK_Nterm')
            else:
                istryptic.append("Error")
            length_end = len(istryptic)
            if length_end == length_begin:
                istryptic.append("Error")
        except:
            istryptic.append("Error")
            continue
    df["istryptic"] = istryptic

    return df


def extend_peptides(df,colname_seq = "UPSEQ",colname_pepseq="PEP.StrippedSequence",size = 3):

    seq = df[colname_seq].to_list()
    pep_seq = df[colname_pepseq].to_list()
    new_pep = []
    for i,x in enumerate(pep_seq):

        m = re.search(str(x),str(seq[i]))
        if m:
            newstart = m.start()-size
            newend = m.end()+size
            if newstart < 0 :
                newstart = 0
            if newend > len(seq[i]):
                newend = len(seq[i])
            new_pep.append(seq[i][newstart:newend])

        else:
            new_pep.append(x)

    df["extended_pep"] =new_pep
    return df

def pep_positioning_along_sequence(df,colname_seq = "UPSEQ",colname_pepseq="PEP.StrippedSequence"):

    seq = df[colname_seq].to_list()
    pep_seq = df[colname_pepseq].to_list()
    relative_pos = []
    for i,x in enumerate(pep_seq):

        m = re.search(str(x),str(seq[i]))
        if m:
            relative_pos.append(m.end()/len(seq[i]))
        else:
            new_pep.append(-1)

    df["relative_position"] =relative_pos
    return df

def compute_ss_peptide(df,colname_ids = "PG.ProteinAccessions",colname_pepseq="PEP.StrippedSequence"):
    with open('/Users/fsesterhenn/Dropbox/PostDoc/databases/psipred/UPID_PSIPRED_human.json', 'r') as fp:
        UPPSI = json.load(fp)
    peplist = df[colname_pepseq].to_list()
    idlist = df[colname_ids].to_list()
    out = []
    for x in idlist:
        try:
            out.append(UPPSI[x])
        except KeyError:
            out.append("NaN")

    df["PSIPRED"] = out



    pep_psi = []
    seqs = df["UPSEQ"].to_list()
    psi = df["PSIPRED"].to_list()
    for i,x in enumerate(peplist):
        m = re.search(str(x),str(seqs[i]))
        if m:
            if str(psi[i]) !="NaN":
                ppsi = str(psi[i][m.start():m.end()])
                pep_psi.append(ppsi)
            else:
                pep_psi.append("X")
        else:
            pep_psi.append("X")



    df["peptide_psipred"] = pep_psi
    df["%loop_pep"] = df["peptide_psipred"].apply(lambda x: x.count("C")/len(x))
    df["%helix_pep"] = df["peptide_psipred"].apply(lambda x: x.count("H")/len(x))
    df["%beta_pep"] =df["peptide_psipred"].apply(lambda x: x.count("E")/len(x))
    return df

def average_replicates_peptide_level(df,colname_pep="PEP.StrippedSequence",colname_intensity="FG.Quantity_y",merge=True):
    ''''this function may need work and needs to be checked properly'''''

    grouped = df.groupby(colname_pep).mean().reset_index()
    grouped_stddev = df.groupby(colname_pep).std().reset_index()
    grouped2 = grouped[[colname_pep,colname_intensity]]
    grouped_stddev2 = grouped_stddev[[colname_pep,colname_intensity]]
    new_colname_intensity = colname_intensity+"_y"
    old_colname_intensity = colname_intensity+"_x"
    if merge:
        df_out = pd.merge(df,grouped2,how='left',on=colname_pep)
        df_out.drop_duplicates(subset=colname_pep,inplace=True)
        df_out.rename({new_colname_intensity:"mean_quantity"},axis=1,inplace=True)
        df_out.drop(old_colname_intensity,axis=1,inplace=True)
        df_out2 = pd.merge(df_out,grouped_stddev2,how='left',on=colname_pep)
        df_out2.drop_duplicates(subset=colname_pep,inplace=True)
        df_out2.rename({"FG.Quantity_y":"stddev_quantity"},axis=1,inplace=True)
        #df_out2.drop(old_colname_intensity,axis=1,inplace=True)

        return df_out2
    else:
        return grouped2

def create_cleavage_region(df,n=False,c=False):

    allaa = ["G","V","L","A",'M',"I",'S','T','C','D','R','E','K','Q','H','N','P',"W","F","Y"]

    peplist = df["extended_pep"].to_list()

    if n:
        for p in peplist:
                try:
                    p2 = p[-6:-1]
                    print(p2)
                except:
                    continue
                    print(p2)
    if c:
        for p in peplist:
                try:
                    p2 = p[0:5]
                    print(p2)
                except:
                    continue
                    print(p2)
