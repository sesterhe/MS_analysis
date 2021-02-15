import os
import numpy
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import *
sys.path.append("/Users/fsesterhenn/local/gits/MS_analysis")
from basic_functions import *
from pdb_functions import *
from database_processing_functions import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.spatial import distance
import re


with open("/Users/fsesterhenn/Dropbox/PostDoc/databases/UNIPROT_functional_annotations/UPID_sequences.json") as fp:
    UPseq  = json.load(fp)


def create_array_foldchange(upids,acc_pep,pep_fc,alignment=False):
    '''''
    this function takes as input a list with uniprot ids. it needs a dictionary with the upids:[peptides] and a dictionary with peptide:foldchange. Faster with alignment=False, will then just perform string matching using regular expression.
    acc_pep = convert_df_to_dict(df3,"PG.ProteinAccessions_x","PEP.StrippedSequence_x")
    pep_acc = convert_df_to_dict(df3,"PEP.StrippedSequence_x","PG.ProteinAccessions_x")
    pep_fc_rapa = convert_df_to_dict(df3,"PEP.StrippedSequence_x","log2_FC_rapa_TC_PK")
    for k,v in pep_fc_rapa.items():
        pep_fc_rapa.update({k:v[0]})

    '''''


    pep_array={}
    for upid in upids:
        pep = acc_pep[upid]
        try:
            ar1 = np.zeros(len(UPseq[upid]))
            for p in pep:
                if alignment:
                    align = align_pep_upseq(p,UPseq[upid])
                    if align[0] >0.9*len(p):
                        start = align[1]
                        end = align[2]
                        ar1[start:end] = pep_fc[p]
                        #print(ar1)
                else:
                    match = re.search(p,UPseq[upid])
                    ar1[match.start():match.end()] = pep_fc[p]


            pep_array.update({upid:ar1})

        except:
            print(upid+ " failed")

    return pep_array


def compute_distance_fingerprints(cond1,cond2):
    output = {}
    for k in cond1.keys():
        dst = distance.euclidean(cond1[k], cond2[k])
        output.update({k:dst})
    return output

def compute_distance_fingerprints2(cond1,cond2):
    output = {}
    for k in cond1.keys():
        if k in cond2.keys():
            a = numpy.nan_to_num(cond1[k],posinf=0, neginf=0)
            b= numpy.nan_to_num(cond2[k],posinf=0, neginf=0)
            dst = distance.euclidean(a, b)
            output.update({k:dst})

    return output

def visualize_fingerprint(array,UPID,cmap='Blues',save=False):

    dpi = 150
    fig = plt.figure(figsize=(4, 2), dpi=dpi)
    plt.imshow(array[UPID].reshape(1, -1),vmin =0, vmax=10 ,cmap=cmap, aspect='auto',interpolation='none')
    plt.xticks()
    plt.rcParams['font.family'] = 'sans-serif'
    plt.xticks(fontsize=14, fontname="Arial")
    plt.tick_params(axis="x",which='major')
    plt.xlabel("Sequence",fontname="Arial",fontsize=14)
    plt.yticks([])
    plt.tight_layout()
    if save:
        outfilename=UPID+"_fingerprint.pdf"
        plt.savefig("/Users/fsesterhenn/Desktop/"+outfilename)
