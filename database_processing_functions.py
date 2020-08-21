

import rstoolbox.io as rs
import json


def generate_dict_from_fasta(infile,outfile):
    '''''
    this function was written to generate a python dictionary from a uniprot fasta file. It will return a dictionary with the UPIDs as keys and the sequence as values, stored in a json file.

    USAGE: generate_dict_from_fasta("/Users/fabian/Downloads/uniprot_sprot.fasta","output.json")

    The json file can later be loaded again by:
    with open('output.json', 'r') as fp:
    data = json.load(fp)
    '''''
    df = rs.read_fasta(infile)
    names = df["description"].to_list()
    newnames = []
    for n in names:
        n=n.split("|"[0])[1]
        newnames.append(n)
    seqs = df["sequence_A"].to_list()
    if len(newnames) == len(seqs):
        dic = dict(zip(newnames, seqs))
        with open(outfile, 'w') as fp:
            json.dump(dic, fp)
    else:
        raise Exception('length do not match')
    return dic



def convert_df_to_dict(df, key, value):
    '''''
    This function takes as input a dataframe and returns a dictionary of two selected columns. The values are organized as a list with possibly multiple members.
    '''''
    dic = {k: list(set(g[value].tolist())) for k,g in df.groupby(key)}
    return dic
