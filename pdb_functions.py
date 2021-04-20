import json
import pandas as pd
import os
import re
import numpy as np
from Bio.Seq import Seq
from Bio import pairwise2
import Bio.PDB


def retrieve_pdb(inputdict):
    import os
    url = "https://files.rcsb.org/download/"
    pdblist = []
    chainlist = []
    namelist =[]

    for e in inputdict.keys():
        try:
            pdb = UPPDB[e][0]

        except KeyError:
            continue

        pdb2 = pdb.split("_")[0]+".pdb"
        chain = pdb.split("_")[1]

        cmd = "wget "+url+pdb2

        if not os.path.isfile(pdb2):
            print(cmd)
            os.system(cmd)

    return None

def retrieve_pdb2(name):
    url = "https://files.rcsb.org/download/"

    cmd = "wget "+url+name

    if not os.path.isfile(name):
        print(cmd)
        os.system(cmd)

def extract_sequence_from_pdb(infile, chain):
    threetoone = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    pdb = open(infile,"r")
    pdb = pdb.readlines()
    sequence = []
    for line in pdb:
        if line.startswith("ATOM"):
            if line[21:22].strip() == chain:
                if line[13:16].strip() == "CA":
                    sequence.append((threetoone[line[17:20]]))

    return("".join(sequence))

def get_indices(seq,peptide):
    start = str(seq).index(peptide)
    length = len(peptide)
    end = start + length
    return [start,end]

def parse_B_factor(infile,chain):
    '''''
    This function is able to get the B factor column from a PDB file. It needs the PDB file downloaded in the same location.
    '''''
    pdb = open(infile,"r")
    pdb = pdb.readlines()
    bfactors = []
    for line in pdb:
        if line.startswith("ATOM"):
            if line[21:22].strip() == chain:
                if line[13:16].strip() == "CA":
                    bfactors.append(float(line[61:66]))
    return bfactors

def retrieve_B_factors_of_selection(dic_msoutput):
    '''''
    This function will retrieve the B factors (according to parse_B_factor function) for a set of peptides.
    The input is a dictionary containing uniprot IDs:peptides.
    The output is a dictionary with the pdb as key and a nested list containing the average B factor of the
    structure followed by the B factor of the peptide as values.


    '''''
    from statistics import mean

    output = {}

    for e in dic_msoutput.keys():
        out = []
        try:
            pdb = UPPDB[e][0].strip()
            #print(pdb)

        except KeyError:
            continue
        if pdb:
            pdb2 = pdb.split("_")[0]+".pdb"
            chain = pdb.split("_")[1]
            s = extract_sequence_from_pdb(pdb2,chain)
            #print(s)

            try:
                #print(dic_msoutput[e][0])

                if dic_msoutput[e][0].strip() in s:
                    ind = get_indices(s,dic_msoutput[e][0].strip())
                    #print(ind)
                    l = parse_B_factor(pdb2,chain)
                    #print(l[ind[0]:ind[1]])
                    out.append(mean(l))
                    out.append(mean(l[ind[0]:ind[1]]))
                    #print(out)
                    #print(pdb2+" mean_B_factor: "+mean(l))
                    #print(dic_msoutput[e][0]+" peptide_B_factor: "+mean(l[ind[0]:ind[1]]))
                    output.update({pdb:out})
                else:
                    print("peptide not found in seq")
            except:
                continue

    return output

def get_ss_of_peptides(query):

    '''''
    input is a dictionary with format UPID:peptide

    '''''
    ss_assignment = {}
    for s in query.keys():
        pdb = []
        try:
            pdb.append(UPPDB[s])
        except KeyError:
            e1 = "no corresponding PDB for "+s
            #print(e1)
            continue
        pdb2 = [val for sublist in pdb for val in sublist]
        if len(pdb2) == 0:
                #print("length 0")
                continue
        else:
            for i,p in enumerate(pdb2):
                ss = []
                #print(pdb2)
                #print(i,p)
                try:
                    p = p.strip()
                    seq = PDBseqs[str(p)]
                except KeyError:
                    #print("failed")
                    continue
                try:
                    p = p.strip()
                    dssp = DSSP[str(p)]
                except KeyError:
                    #print("no dssp available")
                    continue

                for x in query[s]:
                        #print(x)
                        if str(x) in str(seq):
                                start = str(seq).index(x)-1
                                length = len(x)
                                end = start + length
                                if len(ss)>1:
                                    continue
                                else:
                                    ss.append(dssp[start:end])

                                #print(dominant_ss)

                                pdbid = pdb2[i].split("_")[0]
                                pdbchain = pdb2[i].split("_")[1]


                        else:
                            e2 = "not found"
                            #print(e2)

            t = get_dominant_ss(ss)
            ss_assignment.update({x:t})
    E = []
    H = []
    L = []
    NA = []
    for y in ss_assignment.values():
        y2 = "".join(y)
        if y2 == "E":
            E.append(y2)
        if y2 == "H":
            H.append(y2)
        if y2 == "x":
            L.append(y2)
        if y2 == "":
            NA.append(y2)

    #E = sum(value == "E" for value in ss_assignment.values())
    print("E: "+str(len(E)))
    print("H: "+str(len(H)))
    print("L: "+str(len(L)))
    print("N.A.: "+str(len(NA)))
    return(ss_assignment)

def get_PDBs_highlight_peptides2(query,fetchall=False,maxlim=5):

    '''''
    input is a dictionary with format UPID:peptide
    output is a file to copy/paste in pymol, containing the commands to retrieve pdb files,
    group them and highlight the peptides by color.
    This function relies on  the find_seq.py script for finding and highlighting the peptide sequence.

    '''''

    log = open("error2.log","w")
    output = open("pymol_commands2.txt","w")
    output.write("run findseq.py"+"\n")
    for s in query.keys():
        pdb = []
        pdb3 = []
        try:
            pdb.append(UPPDB[s])
        except KeyError:
            e1 = "no corresponding PDB for "+s
            log.write(e1+"\n")
            continue
        pdb2 = [val for sublist in pdb for val in sublist]
        if len(pdb2) == 0:
                print("length 0")
                continue
        else:
            for i,p in enumerate(pdb2[:maxlim]):
                #print(pdb2)
                #print(i,p)
                try:
                    p = p.strip()
                    seq = PDBseqs[str(p)]
                except KeyError:
                    print("failed")
                    continue
                    #sequence_offset = int(("".join(offset[p])))
                command1 = "fetch " + p

                for x in query[s]:
                        #print(x)
                        if str(x) in str(seq):
                                start = str(seq).index(x)-1
                                length = len(x)
                                end = start + length
                                pdbid = pdb2[i].split("_")[0]
                                pdbchain = pdb2[i].split("_")[1]
                                find = "findseq "+x+", "+p+", "+p+"_"+x
                                coloring  = "color red, "+p+"_"+x
                                pdb3.append(p)
                                print(command1+"\n")
                                print(find)
                                print(coloring+"\n")
                                output.write(command1+"\n")
                                output.write(find+"\n"+coloring+"\n")
                                output.write("delete "+p+"_"+x+"\n")
                        else:
                            e2 = "not found"
                            log.write(e2+"\n")
                            if fetchall:
                                output.write(command1+"\n")
        if fetchall:
            output.write("group "+ s+","+(" ".join(x for x in pdb2))+"\n")
            print("group "+ s+","+(" ".join(x for x in pdb2)))
        else:
            output.write("group "+ s+","+(" ".join(x for x in pdb3))+"\n")
            print("group "+ s+","+(" ".join(x for x in pdb3)))
    log.close()
    output.close()
    return log, output

def reduce_ss(ss):
    ss_reduced = []
    for s in ss:
        s =s.replace("G","H")
        s =s.replace("I","H")
        s =s.replace("B","E")
        s =s.replace("T","x")
        s =s.replace("S","x")
        s =s.replace("C","x")
        s =s.replace("L","x")
        ss_reduced.append(s)
    sss = "".join(ss_reduced)
    return sss

def extract_CA_coordinates(infile,chain):
    '''''
    This function is made to get the CA coordinates from a PDB file. It needs the PDB file downloaded in the same location.The output is a file ending on _CA.pdb with the CA coordinates only.
    '''''
    pdb = open(infile,"r")
    pdb = pdb.readlines()
    CA_coord = []
    outfile = infile.strip(".pdb")+"_"+chain+"_CA.pdb"
    out = open(outfile,"w")
    for line in pdb:
        if line.startswith("ATOM"):
            if line[21:22].strip() == chain:
                if line[13:16].strip() == "CA":
                    #CA_coord.append(line)
                    out.write(line)
                    print(line)
    return out


def align_pep_PDBseq(peptide,PDB,chain):
    '''''
    aligns a peptide sequence to the PDBsequence extracted by "extract_sequence_from_pdb" function. Needs peptide,PDB and chain info.
    Returns a list of len 3, where element 1 is an alignment score, element 2 and 3 are the start and end positions of the alignments.
    USAGE: align_pep_PDBseq("RGHIY","2B14.pdb","A")

    '''''
    pdb_seq = extract_sequence_from_pdb(PDB,chain)
    alignments = pairwise2.align.localms(pdb_seq, peptide,1, -1, -1, 0)
    score = alignments[0][2]
    start = alignments[0][3]
    stop = alignments[0][4]
    return [score,start,stop]

def extract_HETATM_coordinates(infile,chain):
    '''''
    This function is able to get the HETATM (except water) coordinates from a PDB file. It needs the PDB file downloaded in the same location. The output is a file ending on _HETATM.pdb with the coordinates.
    '''''
    i = open(infile,"r")
    i = i.readlines()
    outfile = infile.strip(".pdb")+"_HETATM.pdb"
    o = open(outfile,"w")
    for line in i:
        if line.startswith("HETATM"):
            if line[17:20] != "HOH":
                print(line.strip())
                o.write(line)
    #return o.write(line)

def extract_HETATM_ID(infile,chain):
    '''''
    This function is able to get the HETATM (except water) ID from a PDB file. It returns a dictionary with the format PDB:[NAMES OF LIGANDS].

    '''''

    i = open(infile,"r")
    i = i.readlines()
    ligands = []
    ligands2 = []
    PDBLIG={}
    for line in i:
        if line.startswith("HETATM"):
            if line[17:20] != "HOH":
                ligand = line[17:20]
                ligands.append(ligand)
    ligands2 = list(set(ligands))
    PDBLIG[infile.strip("_HETATM.pdb")] = ligands2
    return PDBLIG

def compute_minimal_distance_marker_to_HETATM(marker,HETATM):

    '''''
    computes the minimal distance in A between two pdb files. Input should be a file with the CA coordinates (but can be anything)
    of an MS detected peptide, followed by a coordinate file containing only the HETATMs of a PDB of interest. returns the minimal
    distance.
    '''''

    i = open(marker,"r").readlines()
    coords_marker = []
    distances = []
    coords_hetatm = []
    for x in i:
        #print(x)
        x_coord = float(x[31:38].strip())
        y_coord = float(x[39:46].strip())
        z_coord = float(x[47:54].strip())
        ar = np.array([x_coord,y_coord,z_coord])
        coords_marker.append(ar)


    i = open(HETATM,"r").readlines()
    coords_hetatm = []
    for x in i:
        try:
            x_coord = float(x[31:38].strip())
            y_coord = float(x[39:46].strip())
            z_coord = float(x[47:54].strip())
            ar = np.array([x_coord,y_coord,z_coord])
            coords_hetatm.append(ar)
        except ValueError:
            continue

    for c in coords_marker:
            for c2 in coords_hetatm:
                dist = np.linalg.norm(c-c2)
                distances.append(dist)

    #print(min(distances))
    return(min(distances))

def extract_CA_coordinates_of_marker(msdict):
    '''''
    needs prior loading of UPPDB dictionary. Will load the dictionary with the uniprot IDs and the MS peptides automatically.
    Will generate pdb files with the CA coordinates of the identified peptide(s).
    needs further testing
    USAGE: extract_CA_coordinates_of_marker("/Users/sesterhe/Dropbox/PD_project/PD_vs_healthy/biomarker_dict.json")
    '''''

    with open(msdict, 'r') as fp:
        biomarker = json.load(fp)
    for k in biomarker.keys():
        try:
            #print(k)
            pdb = UPPDB[k][0]
            #print(pdb)
        except KeyError:
            continue
        if pdb:

            pdb_name = pdb.split("_")[0]+".pdb"
            chain = pdb.split("_")[1]
            extract_CA_coordinates(pdb_name,chain)
            seq = extract_sequence_from_pdb(pdb_name,chain)
            #print(seq)
            pep= biomarker[k]
            #print(pep)

            for p in pep:
                try:
                    indices = get_indices(seq,p)
                    i = open(pdb_name.strip(".pdb")+"_"+chain+"_CA.pdb","r")
                    o = open(pdb_name.strip(".pdb")+"_marker_CA.pdb","w")
                    i2 = i.readlines()
                    #print(i2)
                    for e in i2[indices[0]:indices[1]]:
                        #print(e)
                        o.write(e.strip()+"\n")
                except ValueError:
                    continue


#this function needs to be revised to work properly with the numbering
def get_PDBs_highlight_peptides(query):

    '''''
    THIS FUNCTION MAY NEED WORK
    input is a dictionary with format UPID:peptide
    output is a file to copy/paste in pymol, containing the commands to retrieve pdb files,
    group them and highlight the peptides by color.
    This function relies on residue numbering, assuming that the PDB does not start with negative numbers.
    This may cause problems, so be careful to check that the correct sequence is highlighted.
    If according to the info in the pdb_chain_uniprot.csv file there is a discrepany in the
    numbering between PDBSEQ records and numbering of the sequence, it will correct for it by
    renumbering the sequence, but no 100% guarantee that it works in every case.
    An alternative could be to use the find_seq.py script for these cases.

    '''''

    log = open("error.log","w")
    output = open("pymol_commands.txt","w")
    for s in query.keys():
        pdb = []
        try:
            pdb.append(UPPDB[s])
        except KeyError:
            e1 = "no corresponding PDB for "+s
            log.write(e1+"\n")
            continue
        pdb2 = [val for sublist in pdb for val in sublist]
        if len(pdb2) == 0:
                print("length 0")
                continue
        else:
            for i,p in enumerate(pdb2):
                #print(pdb2)
                #print(i,p)
                try:
                    p = p.strip()
                    seq = PDBseqs[str(p)]
                    sequence_offset = int(("".join(offset[p])))
                    command1 = "fetch " + p
                    print(command1+"\n")
                    output.write(command1+"\n")
                    if sequence_offset !=1:
                        command2 = "alter " + p +",resi=str(int(resi)+"+str((int(sequence_offset)-1))+")"
                        print(command2+"\n")
                        output.write(command2+"\n")
                except KeyError:
                    print("failed")
                    continue

                for x in query[s]:
                        #print(x)
                        if str(x) in str(seq):
                                start = str(seq).index(x)-1

                                length = len(x)
                                end = start + length
                                pdbid = pdb2[i].split("_")[0]
                                pdbchain = pdb2[i].split("_")[1]

                                coloring  = "color red, " + p +" and chain "+  pdbchain +" and resi "+str(start)+"-"+str(end)
                                print(coloring+"\n")
                                output.write(coloring+"\n")
                        else:
                            e2 = "not found"
                            log.write(e2+"\n")
        output.write("group "+ s+","+(" ".join(x for x in pdb2))+"\n")
        print("group "+ s+","+(" ".join(x for x in pdb2)))
    log.close()
    output.close()
    return log, output

def find_representative_pdb(pep_dict):
    '''''
    identifies a representative pdb structure or homology model (one, not several). Needs UPPDB_with_homol_models dictionary preloaded. Extracts sequences on-the-fly.
    Input is a pep_dict, generated by: 'pep_dict = convert_df_to_dict(df,"Prot","Pep")'

    returns pep_pdb dictioanary that has the peptide and a representative pdb structure that contains this peptide sequence. can be added to original dataframe
    by df["representative_PDB"] = df["Pep"].map(find_representative_pdb_or_homol(pep_dict))
    '''''

    with open('/Users/fsesterhenn/Dropbox/PostDoc/PD_project/databases/UPID_PDB_CHAIN.json', 'r') as fp:
        UPPDB = json.load(fp)
    pep_pdb_only = {}
    for upid,pep in pep_dict.items():
        count = 0
        try:
            pdb = UPPDB[upid]
        except KeyError:
            continue
        for p in pdb:
            if count < 1:
                p=str(p)
                name = p.split("_")[0]+".pdb"
                chain = p.split("_")[1]
                if not os.path.isfile(name):
                    retrieve_pdb2(name)
                try:
                    ss = extract_sequence_from_pdb(name,chain)
                except FileNotFoundError:
                    continue
                if ss:
                    for pe in pep:
                        if pe in ss:
                            count = 1
                            pep_pdb_only.update({pe:p})

            else:
                continue
    return pep_pdb_only

def find_representative_pdb_or_homol(pep_dict):
    '''''
    identifies a representative pdb structure or homology model (one, not several). Needs UPPDB_with_homol_models dictionary preloaded. Extracts sequences on-the-fly.
    Input is a pep_dict, generated by: 'pep_dict = convert_df_to_dict(df,"Prot","Pep")'

    returns pep_pdb dictioanary that has the peptide and a representative pdb structure that contains this peptide sequence. can be added to original dataframe
    by df["representative_PDB"] = df["Pep"].map(find_representative_pdb_or_homol(pep_dict))
    '''''
    pep_pdb = {}
    for upid,pep in pep_dict.items():
        count = 0
        try:
            pdb = UPPDB_with_homol_models[upid]
        except KeyError:
            continue
        for p in pdb:
            if count < 1:
                p=str(p)

                if len(p) > 7:
                    name = p
                    chain = ""
                else:
                    name = p.split("_")[0]+".pdb"
                    chain = p.split("_")[1]
                    if not os.path.isfile(name):
                        retrieve_pdb2(name)
                try:
                    ss = extract_sequence_from_pdb(name,chain)
                except FileNotFoundError:
                    continue
                if ss:
                    for pe in pep:
                        if pe in ss:
                            count = 1
                            pep_pdb.update({pe:p})

            else:
                continue
    return pep_pdb


def extract_coords2(infile, chain):
    threetoone = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    pdb2 = "5LPE_B"


    pdb = open(infile,"r")
    chain = chain
    pdb = pdb.readlines()
    sequence = []
    coords_peptide = []

    output = {}
    for line in pdb:
        if line.startswith("ATOM"):
            if line[21:22].strip() == chain:
                if line[13:16].strip() == "CA":
                    sequence.append((threetoone[line[17:20]]))
                    xcoords=float(line[30:38].strip())
                    ycoords=float(line[39:46].strip())
                    zcoords=float(line[47:54].strip())
                    ar = np.array([xcoords,ycoords,zcoords])
                    coords_peptide.append(ar)
    seq = "".join(sequence)
    return coords_peptide

def find_closest_active_site(pep_coord_dict):
    '''''
    needs as input a nested dictionary with {UPID:{pep:coords}}
    '''''
    pep_active_site_distances = {}
    for k,v in pep_coord_dict.items():
        try:
            active = active_site_coord[k]
            #print(active)
        except:
            continue
        if str(active) == "nan":
            continue
        else:

            for kk,vv in v.items():
                #print(kk)
                temp_active = []
                for vvv in vv:

                    for a in active:
                        dist = np.linalg.norm(vvv-a)
                        temp_active.append(dist)



                pep_active_site_distances.update({kk:min(temp_active)})
    return pep_active_site_distances

def find_closest_binding_site(pep_coord_dict):
    '''''
    needs as input a nested dictionary with {UPID:{pep:coords}}
    '''''
    pep_binding_site_distances = {}
    for k,v in pep_coord_dict.items():
        try:
            binding = binding_site_coord[k]
            #print(active)
        except:
            continue
        if str(binding) == "nan":
            continue
        else:

            for kk,vv in v.items():
                #print(kk)
                temp_binding = []
                for vvv in vv:

                    for b in binding:
                        dist = np.linalg.norm(vvv-b)
                        temp_binding.append(dist)
                        #print(dist)



                pep_binding_site_distances.update({kk:min(temp_binding)})
    return pep_binding_site_distances
        #print(k,kk, min(temp_binding))
def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain) :
    """Returns a matrix of C-alpha distances between two chains"""
    dist_matrix = numpy.zeros((len(chain), len(chain)), numpy.float)
    for row, residue_one in enumerate(chain) :
        for col, residue_two in enumerate(chain):
            dist_matrix[row, col] = calc_residue_dist(residue_one, residue_two)
    return dist_matrix

def extract_coordinates_peptides(pep_dict,seq=seq,data_representative_pdb=data_representative_pdb):
    '''''
    Input: dictionary with the UniProtID:[peptides].
    Requirements:
    1) load UniProtID:sequence dictionary ("seq").
    2) Create dictionary with peptide:PDB from your dataframe ("data_representative_pdb")
    Output: dictionary containing the peptides with their coordinates as np.array (pep_coord)

    '''''
    pep_coord = {}
    for k, v in pep_dict.items():
        temp_dict = {}
        try:
            for x in v:

                m = re.search(x,seq[k])
                pdb = data_representative_pdb[x]
                name = pdb.split("_")[0]+".pdb"
                chain = pdb.split("_")[1]

                if not os.path.isfile(name):
                    retrieve_pdb2(name)

                s  = extract_sequence_from_pdb(name,chain)


                o = re.search(x,s)
                xx = extract_coords(name,chain)
                pep_coordinates = xx[o.start():o.end()]
                temp_dict.update({x:pep_coordinates})
                pep_coord.update({k:temp_dict})
                #append_value(pep_coord,k,{x:pep_coordinates})

        except KeyError:
            continue
    return pep_coord

def compute_distance_to_func_site(functionalsites,df,pep_coord,colname_ids =UPID):

    '''''
    input:
    functionalsites = dictionary with uniprot ids and a list of residues annotated as functional site (e.g. binding, active site)
    EXAMPLE:
    with open('./data/ORGANIMS/ORGANISM_active_sites.json', 'r') as fp:
        activesites = json.load(fp)
    df = your dataframe
    colname_ids = name of the column that contains the UniProtIDs
    pep_coord = output from extract_coordinates_peptides() function

    output:
    dictionary with peptide:minimal distance to a functionalsite (minimal distance (CAlpha) to the closest functional site only!)

    EXAMPLE USAGE:
    act_dist = compute_distance_to_func_site(activesites,df,p,UPID)
    '''''

    # 1) Create artificial peptides around the functional site (+/-3 residues)

    active_site_peptides = {}
    for k,v in functionalsites.items():
        if str(v) == 'nan':
            continue
        else:
            try:
                seqs = seq[k]
            except KeyError:
                continue
            pept = []
            for x in v:
                pept.append(seqs[x-3:x+3])
                active_site_peptides.update({k:pept})

    # 2) Extract coordinates from PDB for the functional site of interest

    active_site_coord = {}
    UPID_representative_PDB = convert_df_to_dict(df,colname_ids,"representative_PDB")

    for k, v in active_site_peptides.items():

        try:
            representative_pdb = UPID_representative_PDB[k][0]
        except KeyError:
            continue
        if str(representative_pdb) =="nan":
            continue
        else:
            name = representative_pdb.split("_")[0]+".pdb"
            chain = representative_pdb.split("_")[1]
            s  = extract_sequence_from_pdb(name,chain)
            #print(k,s)

            for x in v:
                try:
                    match = re.search(x,s)
                    match2 =  match.start()+2
                    xx = extract_coords(name,chain)
                    pep_coordinates = xx[match2]
                    append_value(active_site_coord,k,pep_coordinates)


                except:
                    continue

    # 3) Compute minimal distance distance

    pep_min_dist_to_active_site = {}
    for k,v in pep_coord.items():
        try:
            active = active_site_coord[k]
        except:
            continue
        if str(active) == "nan":
            continue
        else:
            for kk,vv in v.items():
                temp_active = []
                for vvv in vv:

                    for a in active:
                        dist = np.linalg.norm(vvv-a)
                        temp_active.append(dist)

                pep_min_dist_to_active_site.update({kk:min(temp_active)})

    return pep_min_dist_to_active_site


def extract_HETATM_ID(infile,chain):
    '''''
    This function is able to get the HETATM (except water) ID from a PDB file. It returns a dictionary with the format PDB:[NAMES OF LIGANDS].
    Input: a PDB file and the chain name.
    Example: extract_HETATM_ID("PATH/TO/FILE/6xvf.pdb","A")
    '''''

    i = open(infile,"r")
    i = i.readlines()
    ligands = []
    ligands2 = []
    PDBLIG={}
    for line in i:
        if line.startswith("HETATM"):
            if line[17:20] != "HOH":
                ligand = line[17:20]
                ligands.append(ligand)
    ligands2 = list(set(ligands))
    PDBLIG[infile.strip("_HETATM.pdb")] = ligands2
    return PDBLIG

def compute_distance(PDB1, chain1, pep1, PDB2, chain2, pep2,mode="CA"):

    '''''
    Function to compute the distance between two peptides in either the same or different PDB files.
    Example: compute_distance("/PATH/TO/FILE/4jhw.pdb","CHAIN","PEPTIDEK","/PATH/TO/FILE/4jhw.pdb","CHAIN","PEPTIDEK",mode="CA")
    currently only works in Calpha mode. Measuring distances between all atoms will be added later.

    '''''

    i1 = open(PDB1,"r").readlines()
    i2 = open(PDB2,"r").readlines()
    coords_pep1 = []
    coords_pep2 = []
    distances = []
    i1_seq = extract_sequence_from_pdb(PDB1, chain1)
    i2_seq = extract_sequence_from_pdb(PDB2, chain2)

    i1_list = []
    i2_list = []


    match1 = re.search(pep1,i1_seq)
    if match1 is None:
        raise Exception("peptide sequence "+pep1 + " was not found in "+PDB1)


    match2 = re.search(pep2,i2_seq)
    if match2 is None:
        raise Exception("peptide sequence "+pep2 + " was not found in "+PDB2)


    if mode=="CA":

        for line in i1:
            if line.startswith("ATOM"):
                if line[21:22].strip() == chain1:
                    if line[13:16].strip() == "CA":
                        i1_list.append(line)

        for x in i1_list[match1.start():match1.end()]:

            x_coord = float(x[31:38].strip())
            y_coord = float(x[39:46].strip())
            z_coord = float(x[47:54].strip())
            ar = np.array([x_coord,y_coord,z_coord])
            coords_pep1.append(ar)


        for line in i2:
            if line.startswith("ATOM"):
                if line[21:22].strip() == chain2:
                    if line[13:16].strip() == "CA":
                        i2_list.append(line)

        for x in i2_list[match2.start():match2.end()]:

            x_coord = float(x[31:38].strip())
            y_coord = float(x[39:46].strip())
            z_coord = float(x[47:54].strip())
            ar = np.array([x_coord,y_coord,z_coord])
            coords_pep2.append(ar)

        for c in coords_pep1:
            for c2 in coords_pep2:
                dist = np.linalg.norm(c-c2)
                distances.append(dist)


    print("Minmal distance is: "+str(min(distances)))

    return min(distances)

def compute_distance_HETATM(PDB1, chain1, pep1, PDB2, HETATM):

    '''''
    Function to compute the distance between a peptide and a ligand/HETATM.
    Example: compute_distance("/PATH/TO/FILE/4jhw.pdb","CHAIN","PEPTIDEK","/PATH/TO/FILE/4jhw.pdb","LIGAND_ID")
    You can look up the ID of the ligand by calling the function "extract_HETATM_ID(PDBFILE,CHAIN)".
    '''''

    i1 = open(PDB1,"r").readlines()
    i2 = open(PDB2,"r").readlines()
    coords_pep1 = []
    coords_pep2 = []
    distances = []
    i1_seq = extract_sequence_from_pdb(PDB1, chain1)
    match1 = re.search(pep1,i1_seq)
    i1_list = []
    i2_list = []

    if match1 is None:
        raise Exception("peptide sequence "+pep1 + " was not found in "+PDB1)


    for line in i1:
        if line.startswith("ATOM"):
            if line[21:22].strip() == chain1:
                if line[13:16].strip() == "CA":
                    i1_list.append(line)

    for x in i1_list[match1.start():match1.end()]:

        x_coord = float(x[31:38].strip())
        y_coord = float(x[39:46].strip())
        z_coord = float(x[47:54].strip())
        ar = np.array([x_coord,y_coord,z_coord])
        coords_pep1.append(ar)


    for line in i2:
        if line.startswith("HETATM"):

                if line[17:20].strip() == HETATM:
                        i2_list.append(line)

    for x in i2_list:

        x_coord = float(x[31:38].strip())
        y_coord = float(x[39:46].strip())
        z_coord = float(x[47:54].strip())
        ar = np.array([x_coord,y_coord,z_coord])
        coords_pep2.append(ar)

    for c in coords_pep1:
        for c2 in coords_pep2:
            dist = np.linalg.norm(c-c2)
            distances.append(dist)


    print("Minmal distance to ligand is: " + str(min(distances)))

    return min(distances)


# Identify peptides in spatial proximity to any pair of PDB:list of peptides
def identify_proximal_peptides(pep_pdb,distance=11):

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
            path = "/Volumes/biol_bc_picotti_1/Fabian/PDB_database/"+pathdict[pdb_name]+"/"+pdb_name
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
        except FileNotFoundError:
            continue

        for peptide in pep_pdb[p]:
            match = re.search(peptide,pdbseq)
            contact_pos = np.where(np.any(contact_map[match.start():match.end()], axis=0))
            pep_contacts.update({peptide:contact_pos})
            pep_match.update({peptide:list(range(match.start(),match.end()))})


    for p,v in pep_pdb.items():

        for i,vv in enumerate(v):
            for ii in range(i+1,len(v)):
                tocheck = pep_contacts[v[ii]][0]
                if len(intersection(tocheck,pep_match[vv])) > 1:
                    append_value(out_dict,vv,v[ii])
                    append_value(out_dict,v[ii],vv)

    return out_dict


# PFAM domain analysis
def get_mid_position(df,colname_pep,colname_seq):
    '''''
    Function to match peptides to the uniprot sequence and return the middle position.
    Needed to compute whether peptide is in a given domain or not.
    '''''
    peptides = df[colname_pep].to_list()
    sequences = df[colname_seq].to_list()
    positions = []
    for i, p in enumerate(peptides):
        try:
            match = re.search(p,sequences[i])
            mid = int(match.end()-len(p)/2)
            positions.append(mid)
        except:
            positions.append(np.nan)
    df["middle_peptide_position"] = positions
    return df


def is_in_domain(df,colname_id,colname_pep,colname_middle_peptide_position,UP_PFAM_range):
    pep_domain = {}
    pep_upid = convert_df_to_dict(df,colname_pep,colname_id)
    pep_pos = convert_df_to_dict(df,colname_pep,colname_middle_peptide_position)
    for p in pep_pos.keys():
        upid = pep_upid[p][0]
        try:
            pfams = UP_PFAM_range[upid]
            for pf,pos in pfams.items():
                l = list(UP_PFAM_range[upid][pf])
                try:
                    midpos = int(pep_pos[p][0])
                except ValueError:
                    pep_domain.update({p:np.nan})
                if midpos in range(l[0],l[1]):
                    pep_domain.update({p:pf})
                else:
                    pep_domain.update({p:np.nan})
        except KeyError:
            continue
    return pep_domain
