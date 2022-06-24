"""
        
               _ (.".) _    
              '-'/. .\'-'   
    pyCALF      ( o o )     
                 `"-"`  jgs    


        pyCALF stand for python CALcyanin Finder
        The calcyanin protein will be search and annotated within your input file(s) following three steps:
            1) the glycine triplication specific of calcyanin will be searched using HMM profile.
            2) glycine zipper will be annotated individually using specific HMM models for sequences with a glycine triplication.
            3) the N-ter extremity of sequences with a glyX3 will be annotated using blastp and known n-ter as subject database.
"""


import argparse
import os
import sys


import multiprocessing
import logging

import pandas as pd

from .core import annotcter as cter
from .core import annotnter as nter
from .utils import utils as u

DATASDIR = os.path.join(os.path.dirname(__file__), 'datas')

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        #logging.FileHandler("clustering.log"),
        logging.StreamHandler()
    ]
)

print(sys.prefix)

def get_args():
    parser = argparse.ArgumentParser(
        description="""
        pyCALF\n
                           _(__)_        V
                          '-e e -'__,--.__)
                           (o_o)        ) 
                              \. /___.  |
                              ||| _)/_)/
                             //_(/_(/_(

                        
        pyCALF stand for python CALcyanin Finder
        The calcyanin protein will be search and annotated within your input file(s) following three steps:
            1) the glycine triplication specific of calcyanin will be searched using HMM profile.
            2) glycine zipper will be annotated individually using specific HMM models for sequences with a glycine triplication.
            3) the N-ter extremity of sequences with a glyX3 will be annotated using blastp and known n-ter as subject database.
        """,
        formatter_class= argparse.RawTextHelpFormatter
        )
    parser.add_argument('-i', dest='translated_cds_input', required=True, type=str, #action="extend", nargs="+" , 
                        help='yaml file, fasta file or directory containing cds fasta files')
    parser.add_argument('-e', dest='file_extension', type=str, default="faa.gz" , 
                        help='input files extension')            
    parser.add_argument('-o', dest='res_dir', type=str, required=True,
                        help='output directory')   
    #parser.add_argument('--gff', dest='gff_input', type=str,
    #                    help='yaml file, fasta file or directory containing gff files')
    #parser.add_argument('--fna', dest='genome_input', type=str,
    #                    help='yaml file, fasta file or directory containing genomes files (fna)')                
    parser.add_argument('--log', dest='log', type=str,
                        help='log file')                
    parser.add_argument('--glyx3-hmm', dest='glyx3_phmm', type=str, default= DATASDIR + "/GlyX3.hmm" ,
                        help='path to GlyX3 hmm profile (default: %(default)s)" ')                
    parser.add_argument('--domz', dest='domz', type=int, default=10000,
                        help='sequence space size (default: %(default)s)"')     

    parser.add_argument('--glyx3-coverage', dest='gly3_coverage_threshold', 
                        type=int,default=0.62,
                        help="minimal coverage to be considered as a potential calcyanin (default: %(default)s)" )
    # parser.add_argument('--glyx3-qcovhsp', dest='gly3_qcovhps', type=int,default=0)
    parser.add_argument('--glyx3-evalue', dest='gly3_evalue_threshold', 
                        type=int,default=1e-30,
                        help="hit's evalue threshold (default: %(default)s)")
    parser.add_argument('--glyx3-i-evalue', dest='gly3_i_evalue_threshold', 
                        type=int,default=1,
                        help="domain's i_evalue threshold (default: %(default)s)")

    parser.add_argument('--nterdb', dest='nterdb_fa', 
                        default = DATASDIR + "/nterdb.fasta", 
                        help='path to nterdb fasta file (default: %(default)s)')
    parser.add_argument('--nter-mapping-file', dest='nterdb_tab', 
                        default= DATASDIR + "/nterdb.tsv", 
                        help='path to nterdb mapping file (default: %(default)s)')
    parser.add_argument('--nter-coverage', dest='nter_coverage', 
                        type=int,default=80,
                        help="nter minimal coverage (default: %(default)s)")
    parser.add_argument('--nter-evalue', dest='nter_evalue', 
                        type=int,default=1e-07,
                        help="nter evalue threshold (default: %(default)s)")

    parser.add_argument('--gly1-phmm', dest = 'gly1_phmm', 
                        default = DATASDIR + "/Gly1.hmm", 
                        help='path to GlyZip1 hmm profile (default: %(default)s)')
    parser.add_argument('--gly2-phmm', dest = 'gly2_phmm', 
                        default = DATASDIR + "/Gly2.hmm", 
                        help='path to GlyZip2 hmm profile (default: %(default)s)')
    parser.add_argument('--gly3-phmm', dest = 'gly3_phmm', 
                        default = DATASDIR + "/Gly3.hmm", 
                        help='path to GlyZip3 hmm profile (default: %(default)s)')

    parser.add_argument('--glyzip-i-evalue', dest='glyzip_i_evalue', 
                        type=int,default=3.6e-4,
                        help = "glyzip i-evalue threshold (default: %(default)s)" )
    parser.add_argument('--glyzip-evalue', dest='glyzip_evalue', 
                        type=int,default=1,
                        help="glyzip evalue threshold (default: %(default)s)")

    parser.add_argument('--threads', type=int, default = multiprocessing.cpu_count(),
                        help="(default: %(default)s)")
        
    args = parser.parse_args()
    
    return args


def main():
    print(__doc__)
    args = get_args()
    
    res_dir = os.path.abspath(args.res_dir)
     

    os.makedirs(res_dir + "/fastas" , exist_ok=True)
    os.makedirs(res_dir + "/intermediates" , exist_ok=True)
    
    logging.info("loading input sequences :  %s" % args.translated_cds_input)

    sequences = u.easelfasta(args.translated_cds_input)
    
    logging.info("loading GlyX3 HMM profile")
    
    ghmm = u.easelhmm(args.glyx3_phmm)    
    logging.info("Search GlyX3 in sequence database ... ")
    glyx3hits = cter.findglyx3(
        sequences,
        ghmm, 
        glyx3evalue= args.gly3_evalue_threshold, 
        glyx3ievalue = args.gly3_i_evalue_threshold, 
        cpus = args.threads,
        domZ = args.domz
        )
    logging.info("%i sequence with a glycine triplication found." % len(glyx3hits))
    

    if glyx3hits:
        logging.info("make intermediate file for triplication ..." )
        u.maketable(glyx3hits).to_csv(
            res_dir + "/intermediates/calglyx3.csv",
            index=False,
            header=True,
            sep="\t"
        )     
        logging.info("filter input sequences and keep only glyx3+ sequences..." )
        glyx3seqs = cter.filtersequences(
            sequences,glyx3hits,
            args.gly3_coverage_threshold
            )
        logging.info("done : %i glyx3+ sequences"  % len(glyx3seqs))
        if glyx3seqs:

            logging.info("loading glyzip' specific HMM profiles ... " )
            zhmms = [ u.easelhmm(i)  for i in [args.gly1_phmm , args.gly2_phmm , args.gly3_phmm] ]
            logging.info("done.")
            logging.info("search glyzip in glyx3+ sequences")
            glyziphits = cter.findglyzips(
                sequences = glyx3seqs,
                hmms = zhmms,
                glyzipevalue = args.glyzip_i_evalue, 
                cpus = args.threads,
                domZ = args.domz
            ) 
            logging.info("done.")
            # save fasta and table
            logging.info("save intermediates files ... ")
            u.maketable(glyziphits).to_csv(
                res_dir + "/intermediates/calglyzip.csv",
                index=False,
                header=True,
                sep="\t"
            )

            u.writefasta(
                glyx3seqs, 
                res_dir + "/fastas/glyx3seq.fasta"
            )           
            logging.info("done.")
            # start n-ter annotation

            logging.info("search similar N-ter in %s " % args.nterdb_fa )
            df = nter.blastp(                
                res_dir + "/fastas/glyx3seq.fasta",
                args.nterdb_fa,
                args.nter_evalue
            )
            logging.info("done.")

            nterhits = []

            if df is not None and not df.empty:         
                              
                nterhits = nter.nearest_neighboor(
                    df,
                    DATASDIR + "/nterdb.tsv",
                )                
                logging.info("Similar nter found for %i sequences." % len(nterhits) )
            del df
            
            logging.info("save intermediates files ... ")
            u.maketable(nterhits).to_csv(
                res_dir + "/intermediates/calnter.csv",
                index=False,
                header=True,
                sep="\t"
            )
            logging.info("done.")
            logging.info("make feature table ... ")
            df = u.maketable(glyx3hits + glyziphits +  nterhits)
            df.to_csv(
                res_dir + "/features.csv",
                index=False,
                header=True,
                sep="\t"
            )
            logging.info("done.")

        logging.info("Summarizing calcyanin modular organization to %s" % res_dir + "/summary.csv" )
        
        reliable_cpt = 0
        with open(res_dir + "/summary.csv" , "w") as stream: 
            stream.write("accession;flag;nter;cter;is_trusted\n")
            for seqid , featuredf in df.groupby("seqid"):
                cterom = ",".join(
                    list(
                        featuredf[featuredf.desc == "glyzip"].sort_values("start").domid
                        )
                )
                nterom = "".join(list(featuredf[featuredf.desc == "nter"].domid))
                
                flag = u.summarize(nterom,cterom)
                
                verif =  "checked" if flag == "Calcyanin with known N-ter" else "verification required"
                if verif == "checked":
                    reliable_cpt +=1
                stream.write("{};{};{};{};{}\n".format(
                    seqid, flag , nterom, cterom , verif )
                )
        
        logging.info("done.")
        logging.info("%i reliable calcyanin found." % reliable_cpt)
    else:
        logging.info("No calcyanin found ... ")
    logging.info("end")

if __name__ == "__main__":
    sys.exit(main())
