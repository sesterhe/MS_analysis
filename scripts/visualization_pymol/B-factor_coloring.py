#!/usr/bin/env python
# coding: utf-8

# # Change B-factors in PDB file to color according to other parameters

# ##### This script will change the B-factor column for CA atoms in PDB files that subsequently can be loaded into pymol and the structures can be colored according to the new "B-factor". An example case are either intensities from MS data or correlation scores etc - really anything that is a number and that makes sense to highlight in different colors in a structure. Input is a csv file that contains a column with "representative_PDB" also containing the chain information, the peptide and a protein (UniProt) name. The output is, for each provided PDBfile, a new file with the suffix _colored.pdb. ONLY the chain of interest is in the new file, so if you want other chains also load those separately. 

# In[1]:


import pandas as pd
import sys
import os
sys.path.append("/Users/fabian/local/gits/MS_analysis")
from basic_functions import *
from database_processing_functions import *
from pdb_functions import *
pd.set_option('display.max_rows', 1000)
from Bio.Seq import Seq
from Bio import pairwise2


# In[7]:


def align_pep_PDBseq(peptide,PDB,chain):
    pdb_seq = extract_sequence_from_pdb(PDB,chain)
    alignments = pairwise2.align.localms(pdb_seq, peptide,1, -1, -1, 0)
    score = alignments[0][2]
    start = alignments[0][3]
    stop = alignments[0][4]
    return [score,start,stop]


# In[9]:


df = pd.read_csv("/Users/fabian/Dropbox/PostDoc/PD_project/structural_analysis_hits/pdbs/VarPeps_DeltaSD1_VarLiPH0_FT_CookDistConProt_AddMultiMode_151020_analysed.csv")


# In[10]:


Prot_pep = convert_df_to_dict(df,"Prot2","Pep")


# In[11]:


Pep_PDB = convert_df_to_dict(df,"Pep","representative_PDB")


# In[12]:


Prot_PDB = convert_df_to_dict(df,"Prot2","representative_PDB")


# In[13]:


for k,v in Prot_PDB.items():
    correlation_pep = []
    correlation_list = []
    if str(Prot_PDB[k][0]) == "nan" :
        continue
    else:
        filename = k+"_"+Prot_pep[k][0]+".csv"
        if os.path.isfile(filename):
            correlations = pd.read_csv(filename)
            correlations.drop_duplicates(subset="Pep",inplace=True)
            correlations = correlations.sort_values(by="CorNorm")
            correlation_pep = correlations["Pep"].to_list()
            correlation_list = correlations["CorNorm"].to_list()
            correlation_list = [ '%.3f' % elem for elem in correlation_list]
            
            
            pdb_name = Prot_PDB[k][0].split("_")[0]+".pdb"
            pdb_chain = Prot_PDB[k][0].split("_")[1]
            try:
                pdb_seq = extract_sequence_from_pdb(pdb_name,pdb_chain)
            except:
                print("error for "+pdb_name)
                continue
            
            #arbitrary long "empty" list 
            new_b_factors = [-1.01]*50000
            for i,c in enumerate(correlation_pep):
                acceptable_score = len(c)*0.8
                result = align_pep_PDBseq(c,pdb_name,pdb_chain)
                if result[0] < acceptable_score:
                    continue
                else:
                    s = result[1]
                    e = result[2]
                    for x in range(s,e):
                        new_b_factors[x] = correlation_list[i]

                    fi = open(pdb_name,"r")
                    fi = fi.readlines()

                    output = []
                    count = -1
                    newname = pdb_name.strip(".pdb")+"_colored_CorNorm.pdb"
                    fo = open(newname,"w")
                    for line in fi:
                        if line.startswith("ATOM"):
                            if line[21:22].strip() == pdb_chain:
                                if line[13:16].strip() == "CA":
                                    count = count+1
                                    output.append(line[:60]+" "+str(new_b_factors[count])+"           "+line[73:].strip())
                                else:
                                    output.append(line.strip())
                    for l in output:
                        fo.write(l+"\n")


# In[26]:


get_ipython().run_line_magic('pwd', '')


# In[ ]:




