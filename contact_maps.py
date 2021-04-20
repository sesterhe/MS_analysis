
import pandas as pd
import numpy as np
import re
import json


with open("/Users/fsesterhenn/local/gits/gitlab_picotti_eth/structural_analysis_fabian/data/mouse/mouse_active_sites.json") as fp:
    activesites = json.load(fp)

with open("/Users/fsesterhenn/local/gits/gitlab_picotti_eth/structural_analysis_fabian/data/mouse/mouse_bind_sites.json") as fp:
    bindingsites = json.load(fp)

with open("/Volumes/biol_bc_picotti_1-1/Fabian/PDB_database/pathdict.json") as fp:
            pathdict = json.load(fp)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def identify_proximal_functional_sites_from_contact_maps(upid_pep,functionalsites,distance=11,path="/Volumes/biol_bc_picotti_1-1/Fabian/contact_map_predictions/npy_files/"):

    '''''
    input: dictionary containing UPID:list_of_peptides, generated for example by:
    upid_pep = convert_df_to_dict(df,"UPID","PEP.StrippedSequence")

    Make sure you are connected to NAS where most npy files are stored.

    Used are contact maps (reduced distance map with boolean values according to distance cutoff).
    For each peptide in upid_pep, the function will check if a functional site is in proximity (distance cutoff).


    '''''

    with open("/Users/fsesterhenn/Dropbox/PostDoc/databases/UNIPROT_functional_annotations/UPID_sequences.json") as fp:
        UPseq  = json.load(fp)


    pep_proximal_to_func_site = {}

    for p in upid_pep.keys():

        try:
            activesite = functionalsites[p]
            seq = UPseq[p]
            #print(activesite)
            if str(activesite) == "nan":
                continue

        except:
            continue

        try:
            path_file = path+str(p)+".npy"
            print(path_file)
            npz = np.load(path_file)
        except:
            continue


        peplist = upid_pep[p]

        for x in peplist:
            m = re.search(x,seq)
            matchlist = list(range(m.start(),m.end()))

            for a in activesite:

                    contact_map = npz < distance
                    contact_pos = np.where(np.any(contact_map[a-2:a+2], axis=0))[0]

                    for c in contact_pos:

                        if c in matchlist :
                            pep_proximal_to_func_site.update({x:True})
    return pep_proximal_to_func_site


def identify_proximal_peptides_contact_map(upid_pep,distance=11,path="/Volumes/biol_bc_picotti_1-2/Fabian/contact_map_predictions/npy_files/"):

    '''''
    input: dictionary containing UPID:list_of_peptides, generated for example by:
    upid_pep = convert_df_to_dict(df,"upid","PEP.StrippedSequence")

    Make sure you are connected to NAS where npy files are stored.

    The basis for this function are contact maps (reduced distance map with boolean values according to distance cutoff).
    For each peptide in upid_pep, the function will check if other peptides are in proximity (distance cutoff).
    The output is a dictionary with peptides as keys and and the values are all peptides that are in spatial proximity.

    The output can be used to check for each peptide (key) whether any spatially proximal peptide (values) have a fold change or are significant, in order to then mark them as hotspot.

    '''''

    out_dict = {}

    pep_proximal_to_other_signficant_pep = {}

    pep_match = {}

    for upid in upid_pep.keys():

        try:
            path_file = path+str(upid)+".npy"
            npz = np.load(path_file)
        except:
            print("no npy file for "+str(upid))
            continue

        seq = UPseq[upid]

        contact_map = npz < distance

        peplist = upid_pep[upid]


        for pep in peplist:
            contact_positions = []
            m = re.search(pep,seq)
            matchlist = list(range(m.start(),m.end()))
            pep_match.update({pep:matchlist})
            for ma in matchlist:
                contact_pos = np.where(np.any(contact_map[ma-2:ma+2], axis=0))[0].tolist()
                contact_positions.append(contact_pos)
            flatten = [item for sublist in contact_positions for item in sublist]
            contact_positions2 = sorted(list(set(flatten)))
            pep_proximal_to_other_signficant_pep.update({pep:contact_positions2})
            #print(pep_proximal_to_other_signficant_pep)

    for p,v in upid_pep.items():
        for i,vv in enumerate(v):
            try:
                for ii in range(i+1,len(v)):
                    tocheck = pep_proximal_to_other_signficant_pep[v[ii]]
                    if len(intersection(tocheck,pep_match[vv])) > 1:
                        append_value(out_dict,vv,v[ii])
                        append_value(out_dict,v[ii],vv)
            except:
                #print("could not find proximal peptide for "+str(p))
                continue


    return out_dict


def hotspot_spatial_proximity(pep_proximal,pep_FC,cutoff1,cutoff2):

    '''''
    This function will take the output from 'identify_proximal_peptides_contact_map', which is a dictionary
    containing peptides as keys and all spatially proximal peptides as values, and looks up for each peptide key
    if it has a certain foldchange/significance/pivalue (to be adjusted). If the condition is met, it will then
    loop over each value and check whether a second condition (e.g. FC>1.5 or other less stringent criteria) are met.
    pep_FC is a dictionary with the peptide as keys and foldchange or pi-value as value.

    Output: If two or more peptides are spatially proximal and both change (with more or less stringent cutoffs),
    this region will be deemed a Hotspot. The output is a dictionary with the peptide as key and 'True' if its a hotspot.

    USAGE example: df_nostruc["has_proximal_pep_withFC>1.5"] = df_nostruc["pep_stripped_sequence"].map(hotspot_spatial_proximity(pep_proximal,pep_FC,2,1.5))


    '''''


    ishotspot = {}
    for k,v in pep_proximal.items():
        if pep_FC[k][0]>cutoff1:
            if type(v) is list:
                for vv in v:
                    if abs(pep_FC[vv][0])>cutoff2:
                        ishotspot.update({k:True})
            else:
                if abs(pep_FC[v][0])>cutoff2:
                    ishotspot.update({k:True})

    return ishotspot


def identify_proximal_peptides_from_pdb(pep_pdb,distance=11):

    '''''
    input: dictionary containing PDB:list_of_peptides, generated for example by:
    pep_pdb = convert_df_to_dict(df,"representative_PDB","PEP.StrippedSequence")

    Make sure you are connected to NAS where most pdb files are stored.

    Computed are contact maps (reduced distance map with boolean values according to distance cutoff).
    For each peptide in pep_pdb, the function will check if other peptides are in proximity (distance cutoff).
    The output is a dictionary with peptides as keys and and the values are all peptides that are in spatial proximity.

    The output can be used to check for each peptide (key) whether any spatially proximal peptide (values) have a fold change or are significant, in order to then mark them as hotspot.

    '''''


    out_dict = {}
    pep_contacts = {}
    pep_match = {}
    pdb_pep_contacts = {}
    for p in pep_pdb.keys():

        pdb_name = p.split("_")[0]+".pdb"
        pdb_chain = p.split("_")[1]


        try:
            path = "/Volumes/biol_bc_picotti_1-1/Fabian/PDB_database/"+pathdict[pdb_name]+"/"+pdb_name
            fi = open(path,"r")
            fi = fi.readlines()
        except KeyError:
            cmd = "wget https://files.rcsb.org/download/"+pdb_name
            os.system(cmd)
            if not os.path.isfile(pdb_name):
                print("this file was not found")
                continue

        try:
            p2name = p+"_clean.pdb"
            if not os.path.isfile(p2name):
                px = open(path).readlines()
                p2 = open(p2name,"w")
                for l in px:
                    if l.startswith("ATOM"):
                        p2.write(l)
                p2.close()


            structure = Bio.PDB.PDBParser().get_structure(p2name, p2name)
            model = structure[0]
            dist_matrix = calc_dist_matrix(model[pdb_chain])
            contact_map = dist_matrix < distance


            pdbseq = extract_sequence_from_pdb(p2name,pdb_chain)
        except:
            continue

        for peptide in pep_pdb[p]:
            match = re.search(peptide,pdbseq)
            if match:
                contact_pos = np.where(np.any(contact_map[match.start():match.end()], axis=0))
                pep_contacts.update({peptide:contact_pos})
                pep_match.update({peptide:list(range(match.start(),match.end()))})


    for p,v in pep_pdb.items():
            for i,vv in enumerate(v):
                try:
                    for ii in range(i+1,len(v)):
                        tocheck = pep_contacts[v[ii]][0]
                        if len(intersection(tocheck,pep_match[vv])) > 1:
                            append_value(out_dict,vv,v[ii])
                            append_value(out_dict,v[ii],vv)
                except:
                    continue

    return out_dict
