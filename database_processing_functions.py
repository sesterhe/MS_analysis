

#import rstoolbox.io as rs
import json
import statistics
import numpy as np

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


def parse_diso_file(infile):
    av_dis = {}
    mean_dis = {}
    fi  = open(infile).readlines()
    name = infile.split("/")[-1].strip("*.diso")
    disorderscore = []
    symbols = []
    for line in fi:
        if not line.startswith("#"):
            symbols.append(line[8:9])
            disorderscore.append(float(line[10:18].strip()))
    av_disordered_residues = symbols.count("*")/len(symbols)
    m = statistics.mean(disorderscore)
    av_dis.update({name:av_disordered_residues})
    mean_dis.update({name:m})
    with open("UPID_fraction_disorder.json", 'w') as fp:
            json.dump(av_dis, fp)
    with open("UPID_mean_disorder.json", 'w') as fo:
            json.dump(mean_dis, fo)
    return av_dis,mean_dis


def extract_func_site(df,header,pattern,outfile):
    '''''
    extracting functional sites (active site, binding site etc from Uniprot tsv file
    df = pandas dataframe with the info
    header = header of the column of interest (Active site or Binding site or Calcium binding)
    pattern = pattern to search for. Should be ACT_SITE or BINDING for example
    outfile = name of json file to store dictionary with output
    '''''

    l = df[header].to_list()

    entry = df["Entry"].to_list()

    search = pattern

    out = {}

    for i,e in enumerate(l):
        if search in str(e):
            temp = []
            m = re.finditer(search,e)
            for x in m:
                try:
                    temp.append(int(e[x.end():x.end()+5].strip(";").strip()))
                except ValueError:
                    try:
                        temp.append(int(e[x.end():x.end()+4].strip(";").strip()))
                    except ValueError:
                        temp.append(int(e[x.end():x.end()+3].strip(";").strip()))
            print(temp)
            out.update({entry[i]:temp})
        else:
            out.update({entry[i]:e})
    with open(outfile, 'w') as fp:
            json.dump(out, fp)
    return out

def compute_disorder_for_peptides(ids,pep):
    '''''
    this function takes as input a list of uniprot ids and a list of peptides for which disorder should be computed.
    the output is a dictionary with the peptide and the mean disorder score from the precomputed disopred scores, which are loaded from an external hard drive or the nas.
    '''''
    pep_disorder = {}
    for idx,i in enumerate(ids):
        seq = []
        scores = []
        try:
            fi = open("/Volumes/Fabian_data/databases/disorder_all_human/diso_files/"+i+".diso","r")
            fi2 = fi.readlines()
            for l in fi2:
                if not l.startswith("#"):
                    seq.append(l[6:7])
                    #print(l)
                    scores.append(l[10:15].strip())

            sequence = "".join(seq)
            match = re.search(pep[idx],sequence)
            match_scores = scores[match.start():match.end()]
            temp = []
            for m in match_scores:
                temp.append(float(m))
            score = sum(temp)/len(temp)
            pep_disorder.update({pep[idx]:score})
        except FileNotFoundError:
            print("file not found")
            pep_disorder.update({pep[idx]:np.nan})
        #print(pep_disorder)
    return(pep_disorder)
