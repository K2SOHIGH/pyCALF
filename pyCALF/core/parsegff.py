import re
import gffutils

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def load_gff(fn):
    db = gffutils.create_db(fn, dbfn=fn+'.db', force=True, keep_order=True,
        merge_strategy='merge', sort_attribute_values=True)
    return gffutils.FeatureDB(fn+'.db', keep_order=True)

def _feature_id(f):
    if re.search('Prodigal',f.source):
        seqid = f.seqid 
        fid = f.id
        return seqid+"_"+fid.split('_')[-1]
    else:
        fid = "-".join(f.id.split('-')[1:])
        return fid

def gff2dict(gff,feature_types='CDS'):
    gff_dictionnary = {}
    if feature_types:
        genomic_region = ""
        for feature in gff.features_of_type(feature_types):
            if re.search("Prodigal",feature.source):
                if feature.seqid != genomic_region:
                    genomic_region = feature.seqid
                    cpt = 1
                fid = feature.seqid+"_"+str(cpt) 
                cpt+=1            
            else:
                fid = _feature_id(feature)
            if fid not in gff_dictionnary:
                gff_dictionnary[fid]={}
            attributes = []
            for attr in feature.attributes:
                attributes.append(attr + ":" + ";".join(feature.attributes[attr]))        
            gff_dictionnary[fid].update({
                "accession" : feature.seqid,
                "type" : feature.featuretype,
                "start" : feature.start,
                "end" : feature.end,
                "strand" : feature.strand,
                "frame" : feature.frame,
                "src" : feature.source,
                "attr": ",".join(attributes)
            })
    return gff_dictionnary

def check_start_codon(str_seq):
    if str(str_seq)[0:3] in 'ATG':
        return True
    return False

def check_stop_codon(str_seq):
    if str(str_seq)[-3:] in ['TAG','TAA','TGA']:
        return True
    return False

def fasta2dict(fn):
    return SeqIO.to_dict(SeqIO.parse(open(fn),format='fasta'))

def extract_cds_fna(genomic_accession,start,stop,frame,fasta_fna_file):
    fasta = fasta2dict(fasta_fna_file)
    seq = fasta[genomic_accession]
    seq = seq.seq
    if frame == "+":
        start = start-1
        stop = stop
    elif frame == "-":
        start = start - 1
        stop = stop
    else:
        raise ValueError  
    subseq = Seq(seq[start:stop])
    if frame == "-":
        subseq = subseq.reverse_complement()
    if not check_start_codon(subseq):
        print("WARNING NO START CODON DETECTED")
    if not check_stop_codon(subseq):
        print("WARNING NO STOP CODON DETECTED")
    return subseq

def resolve_product_id(accession):
    protid = accession.split('_prot_')[-1].split("_")
    protid.pop(-1)
    return "_".join(protid)
