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
    This function is able to get the CA coordinates from a PDB file. It needs the PDB file downloaded in the same location.
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


    return out


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
