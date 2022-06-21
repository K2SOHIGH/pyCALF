import pandas as pd
from ..utils import utils


# sequences = utils.easelfasta("/Users/maxime/Documents/SRC/modules/pyCALF/pyCALF/datas/GlyX3.msa.fa")


def findglyx3(sequences,hmm, glyx3evalue=1e-30 , glyx3ievalue = 1e-10 ,cpus=4,title="-",**kwargs):
    """pyhmmsearch glyX3 hmm profile in sequences datatabes
    
    Parameter
    ---------
    sequences : list of pyhmmer.easel.DigitalSequence
    hmm : pyhmmer.plan7.HMM
    glyx3evalue: float
        i_evalue threshold 
    cpus : int
    title : str
        extra information regarding set of sequences (e.g assembly accession)
    **kwargs
        will be pass to pyhmmer.hmmsearch pipeline, see documentation for details
                          
    Return
    ------
    list_hit : list
        list of hit as tuple(seqid,hmmname,ndom,evalue,coverage,start,end,title )
        
    """    
    
    #hmmsearch
    glyx3hits = utils.pyhmmsearch(sequences,[hmm],cpus,**kwargs)
    perseqhits = utils.targetview(glyx3hits)    
    list_hit = []

    for seqid, tophits in perseqhits.items():      
        # Only one hit in topHits because only one HMM (i.e GlyX3)        
        if tophits[0].evalue < glyx3evalue:
            start, end = utils.targethitpos(tophits[0] , glyx3ievalue)
        else:
            start, end = utils.targethitpos(tophits[0] , 1000000)

        h = {
            "seqid" : seqid,
            "domid" : hmm.name,
            "start" : start,
            "end" : end,
            "evalue" : tophits[0].evalue,
            "coverage" :  utils.hitcoverage( hmm.M , tophits[0] ),
            "desc": "cter",
            "src": hmm.name,

        }
        
        list_hit.append(
            utils.Hit(**h)            
        )    
    return list_hit


def filtersequences(sequences, datas , coverage_threshold=0.62):
    """filter sequences base on their coverage 
    
    Parameter
    ---------
    sequences : list
        pyhmmer.easel.DigitalSequences    
    datas : list
        list of hit as tuple(seqid,hmmname,ndom,evalue,coverage,start,end,title )
    coverage_threshold: float        
                          
    Return
    ------
    filtered : list
        filtered pyhmmer.easel.DigitalSequences            
    """    
    filtered = []
    
    for i in datas:    
        if i.coverage > coverage_threshold:
            filtered.append(
                 utils.getseqbyname(i.seqid,sequences)
            )
    return filtered



def keep_lowest_evalue(p,hits):
    """keep lowest evalue annoation at a specific position    
    Parameter
    ---------
    p : int
        position within sequence
    hits : list
        sorted list of hits as tuples (evalue , hmmname, ali_from, ali_to , hmm_len) for  the same target (i.e sequence)
    Return
    ------
    dom : str
        hmmname
    evalue: float
        evalue of the associated hit

    """
    # here i do not check the e-value because hits are already sorted.
    # thus if position p is in range of hits h we can break the loop. 
    #h = [seq["i-Evalue"],seq["domains"],seq["alifrom"],seq["alito"]]
    dom = None
    evalue = None
    for h in hits:        
        if p in range(h[2],h[3]):
            dom = h[1]
            evalue = h[0] 
            break #break loop because best hit have been found for p
    return dom,evalue


def filterdomains(domains,threshold=0.65):
    """ remove domain with a coverage above threshold
    Parameter
    ---------
    seqid : str        
    domains : list
        list of domains as tuple
    hmm_len: int
    Return
    ------
    fdomains: list
        list of domains with a coverage above threshold as tuple 
"""    
    fdomains = domains.copy()    
    for i in domains:        
        if i.domid:             
            if i.coverage < threshold:
                fdomains.remove(i)
        else: #linker - no domain 
            fdomains.remove(i)
    return fdomains


def deoverlap(seqid,region,hits):
    """for the same target with multiple domains, for each residue, keep the best annotation with the lowest evalue
    Parameter
    ---------
    seqid : str        
    region : list
        range from sequence start to sequence end 
    hits: list
        sorted list of hits as tuples (evalue , hmmname, ali_from, ali_to , hmm_len) for  the same target (i.e sequence)
    Return
    ------
        list of non-overlaping domains with a coverage above threshold as tuple 

    """
    hmm_len = {h[1]:h[4] for h in hits}
    l_dom = []
    l_evalue = []
    dc = []
    for p in region:
        dom,evalue = keep_lowest_evalue(p,hits)       
        if not l_dom :              # first domain
            l_dom.append(dom)
            l_evalue.append(evalue)
            start=p
        elif dom != l_dom[-1]:      #new domain
            end = p-1               #end of previous domain
            dc.append([start,end])  #append previous born
            l_dom.append(dom)  
            l_evalue.append(evalue)               
            start = p
        else:                       #still the same domain            
            pass    
    end = p
    dc.append([start,end])
    coverage=[]
    for i,j in zip(l_dom,dc):        
        if i:
            coverage.append((j[1]-j[0]) / hmm_len.get(i))
        else:
            coverage.append(None)
    
    # list of non-overlapping domains as tuple (hmmname, [start,end] , evalue, coverage)
    doms = []
    for i in list(zip(l_dom,dc,l_evalue,coverage)):
        doms.append(
                utils.Hit(
                    seqid = seqid,
                    domid = i[0],
                    start = i[1][0],
                    end  = i[1][1],
                    evalue = i[2],
                    coverage = i[3],
                    desc  = "glyzip",
                    src = i[0],
            )
        )

    return filterdomains(doms)


def findglyzips(sequences,hmms,glyzipevalue = 3.6e-4 ,cpus=4,**kwargs):
    """
    Parameter
    ---------
    seqid : str        
    region : list
        range from sequence start to sequence end 
    hits: list
        sorted list of hits as tuples (evalue , hmmname, ali_from, ali_to , hmm_len) for  the same target (i.e sequence)
    Return
    ------
        list of non-overlaping domains with a coverage above threshold as tuple 

    """    
    hmm_length_dict = {i.name:i.M for i in hmms }
    #hmmsearch
    hits = utils.pyhmmsearch(sequences,hmms,cpus,**kwargs)
    
    # del hmm 
    perseqhits = utils.targetview(hits)    

    datas = []
    for seqid, tophits in perseqhits.items():
        seq = utils.getseqbyname(seqid,sequences)        
        pos = range( 0 , len(seq) )
        i_evalue = []
        hmmname = []
        ali_from = []
        ali_to = []
        query_len = [] 

        for hit in tophits:
            for dom in hit.domains :
                if dom.i_evalue < glyzipevalue:                    
                    i_evalue.append(dom.i_evalue)
                    hmmname.append(dom.alignment.hmm_name)
                    ali_from.append(dom.alignment.target_from)
                    ali_to.append(dom.alignment.target_to)
                    query_len.append(hmm_length_dict[dom.alignment.hmm_name])
        h = sorted(list(zip(i_evalue , hmmname, ali_from, ali_to, query_len)))                
        datas += deoverlap(seqid,pos,h)
    return datas