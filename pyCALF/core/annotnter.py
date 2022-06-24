import os
import logging
import subprocess
import pandas as pd
import shutil
import sys
from ..utils import utils as u

def run(cmd):
    """Run the given command line string via subprocess."""
    # Avoid using shell=True when we call subprocess to ensure if the Python
    # script is killed, so too is the child process.
    try:
        child = subprocess.Popen(
            cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except Exception as err:
        sys.exit("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
    stdout, stderr = child.communicate()
    return_code = child.returncode

    # keep stdout minimal as shown prominently in Galaxy
    # Record it in case a silent error needs diagnosis
    if stdout:
        sys.stderr.write("Standard out:\n%s\n\n" % stdout)
    if stderr:
        sys.stderr.write("Standard error:\n%s\n\n" % stderr)

    error_msg = None
    if return_code:
        cmd_str = " ".join(cmd)
        error_msg = "Return code %i from command:\n%s" % (return_code, cmd_str)
    elif "Database or network connection (timeout) error" in stdout + stderr:
        error_msg = "Database or network connection (timeout) error"
    elif "Annotation of 0 seqs with 0 annots finished." in stdout + stderr:
        error_msg = "No sequences processed!"

    if error_msg:
        print(error_msg)
        sys.exit(error_msg)


def blastp(query,subject,evalue = 1e-4):
    import sys
    if sys.version_info[0] < 3: 
        from StringIO import StringIO
    else:
        from io import StringIO


    command = [
            "blastp",
            "-outfmt" , 
            "6 delim=; qacc sacc qlen slen evalue bitscore score pident nident mismatch qstart qend sstart send length qseq sseq qcovs qcovhsp", 
            "-query" , query,
            "-subject", subject,
            "-evalue" , str(evalue),
        ]
        
    o = subprocess.run(command,capture_output=True)
    print(o)
    res = o.stdout.decode('ascii').strip()
    print(res)
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

    