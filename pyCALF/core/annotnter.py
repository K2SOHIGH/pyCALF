import os
import logging
import subprocess
import pandas as pd

from ..utils import utils as u

def blastp(query,subject,evalue = 1e-4):
    import sys
    if sys.version_info[0] < 3: 
        from StringIO import StringIO
    else:
        from io import StringIO
    o = subprocess.run([
            "blastp",
            "-outfmt" , 
            "6 delim=; qacc sacc qlen slen evalue bitscore score pident nident mismatch qstart qend sstart send length qseq sseq qcovs qcovhsp", 
            "-query" , query,
            "-subject", subject,
            "-evalue" , str(evalue),
        ], capture_output=True)
    
    res = o.stdout.decode('ascii').strip()
    if res:
        df = pd.read_table( StringIO(res)  , sep=";"  , header=None )          
        df.columns = "qacc sacc qlen slen evalue bitscore score pident nident mismatch qstart qend sstart send length qseq sseq qcovs qcovhsp".split(" ")
        df["coverage"] = df.apply(lambda x: (x.send-x.sstart) / x.slen * 100, axis=1 )
        return df
    logging.warning("N-ter annotation fail ... ")  
    return None



def load_nter_mapping_file(f:str):    
    nter_db_dict={}
    assert os.path.exists(f)    
    with open(str(f),"r") as db:
        for r in db.readlines():
            nter,seqid = r.strip().split()
            nter_db_dict.update({seqid:nter})
    return nter_db_dict


 
def nearest_neighboor( df:pd.DataFrame , db:str , coverage_threshold:float = 80 , evalue_threshold:float =3.6e-4 ):
    kn = load_nter_mapping_file(db)
    df["domain"] = df["sacc"].map(kn)         

    df = df[(df.evalue <= evalue_threshold) & (df.coverage > coverage_threshold)]                     
    df = df.sort_values(by="evalue",axis=0)
    
    df.drop_duplicates(subset="qacc",keep='first',inplace=True)
    
    nterhits = [] 
    for i,j in df.iterrows():
        nterhits.append(
            u.Hit(
                seqid = j.qacc,
                domid = j.domain, 
                start = j.qstart,
                end = j.qend,
                evalue = j.evalue,
                coverage = j.coverage,
                desc = "nter",
                src = j.sacc,
            ) 
        )
    return nterhits

    