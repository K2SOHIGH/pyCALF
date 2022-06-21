from multiprocessing.sharedctypes import Value
import os
from numpy import select
import pandas as pd

import pyhmmer.easel

class Hit:
    def __init__(self,
        seqid="",
        domid="",
        start:int=0,
        end:int=0,
        evalue:float=0,
        coverage:float=0,
        **kwargs
        ):
        
        
        
        self.seqid = seqid
        self.domid = domid
        self.start = start
        self.end = end
        self.evalue = evalue
        self.coverage = coverage

        if kwargs:
            self.__dict__.update(**kwargs)
            


    def to_dict(self):
        d = {}
        for i,j in vars(self).items():
            if isinstance(j, (bytes, bytearray)):
                d[i] = j.decode('ascii')
            else:
                d[i] = j

        return d


def easelhmm(f:str):
    with pyhmmer.plan7.HMMFile(f) as hmm_file:
        hmm = hmm_file.read()
    return hmm


def easelfasta(f:str,digital:bool=True):
    try:
        with pyhmmer.easel.SequenceFile(f, digital=digital) as seqs_file:
            return list(seqs_file)
    except ValueError:
        raise ValueError("Illegal character found in fasta(s) sequence(s)")


def getseqbyname(name, sequences:list):
    for i in sequences:
        if i.name == name:
            return i
    raise ValueError("Sequence not found : %s" % name)


def pyhmmsearch(sequences:list,hmms:list,cpus=4,**kwargs):
    assert isinstance(hmms,list)
    assert isinstance(sequences,list)    
    
    for hmm in hmms:
        assert isinstance(hmm, pyhmmer.plan7.HMM)
    for seq in sequences:
        assert isinstance(seq, pyhmmer.easel.DigitalSequence)

    all_hits = list(pyhmmer.hmmsearch(hmms,sequences, cpus=4 , **kwargs))
    return all_hits

# def domBeloweEvalue(domains:pyhmmer.plan7.Hits , i_evalue:float = 1e-4):
#     doms = 
#     for d in domains:
#         if d.i_evalue < i_evalue:
#             doms.append(d)
#     return doms


def targetview(tophits):
    """
        hits : list of query hmm tophits
    """
    targets = {}
    for queryhmmhits in tophits:
        for hit in queryhmmhits:
            if hit.name not in targets:
                targets[hit.name] = []
            targets[hit.name].append(hit)
    return targets


def hitcoverage(
    hmm_query_length:int,  
    hit:pyhmmer.plan7.Hit
    ):
    """
        per query (hmm in case of hmmsearch) cumulative coverage.
    """    
    poscovered = []
    for domain in hit.domains:
        poscovered += list(
            range(
                domain.alignment.hmm_from,
                domain.alignment.hmm_to,
            )
        )
    cov = len(list(set(poscovered)))/hmm_query_length
    return cov


def filtertophit(tophits,evalue):
    pass

def targethitpos(hit , evalue = 1e-10):    
    fstart = []
    fend = []          
    
    for d in hit.domains:       
        if d.i_evalue < evalue:
            fstart.append(d.alignment.target_from)
            fend.append(d.alignment.target_to)    
    
    return min(fstart) , max(fend)
    
def _reload():
    import importlib
    importlib.reload(annotcter)
    importlib.reload(annotnter)
    importlib.reload(utils)


def writefasta(sequences:list,file:str):
    f = os.path.abspath(file)
    if not os.path.isdir(os.path.dirname(f)):
        os.makedirs(os.path.dirname(f) , exist_ok = True)
    with open(f,'w') as stream:
        for seq in sequences:
            assert isinstance(seq,pyhmmer.easel.DigitalSequence)
            stream.write(">{} {}\n".format(seq.name.decode('ascii') , seq.description.decode('ascii') ))
            stream.write(seq.textize().sequence + "\n")



def maketable(data:list):
    d = [i.to_dict() for i in data]
    df = pd.DataFrame(d)
    return df

        
def summarize(nter , cter):
    import re
    if re.search("Gly1,Gly2,Gly3",cter):
        if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
            flag = "Calcyanin with known N-ter"
        else:
            flag = "Calcyanin with new N-ter"
    elif re.search("Gly1,Gly3",cter) and nter == "Y-type":
        flag = "Calcyanin with known N-ter"
    elif re.search ("Gly(1|2|3)",cter):
        if nter in ["CoBaHMA-type","Y-type","X-type","Z-type"]:
            flag = "Atypical Gly region with known N-ter"
        else:
            flag = "Atypical Gly region with new N-ter"
    else:
        flag="Ancestral gly containing protein"
    return flag
