# CalcyaninFinder

Find Calcyanin within fasta(s) file(s) in a three steps process :
 
  1) Detecte sequence with a GlyZIp triplication using a HMM profile. 
  2) Resolve C-ter triplication composition using specifics GlyZips profiles. 
  3) Annotate N-ter regions using known N-ter family (X,Y,Z and CoBaHMA) 

## Decision tree : Calcyanin modular organisation and Calcyanin FLAG :

<img src="./resources/img/decision_tree.svg">


## CalcyaninFinder WORKFLOW : 

<img src="./resources/img/calcyanin_finder_workflow.svg">




## Dependencies :
  - Biopython
  - HMMER
  - blast
  - gffutils


## Input :

input_file: yaml file with file identifier as key and file path as value (translated CDS) 
related_gff : yaml file with file identifier as key and file path as value (related GFF file with cds genomic position) 
genome : yaml file with file identifier as key and file path as value (genome nucleotidic sequence(s)) 


## Output :

```
.
├── fastas/ [nter, cter , calcyanin fasta files]
├── GlyX3Finder/
│   ├── assemblies/
│   ├── batch-0
│   │   ├── dom_tables
│   │   └── tables
│   ├── fastas
│   └── tables
├── GlyX3Solver
│   └── tables
├── NterSolver
│   └── tables
└── tables [calcyanin_feature_table.tsv  calcyanin.json  calcyanin.pkl  calcyanin.tsv
```

  - calcyanin_feature_table.tsv  [ feature table ~gff with start and stop position for N-ter and glycine zippers ]
  - calcyanin.json ...
  - calcyanin.pkl  ...
  - calcyanin.tsv ...

json entry example:
```json
{"lcl|NZ_CP089094.1_prot_WP_084990073.1_261": 
  {"partial": "00", 
  "pseudo": false, 
  "assembly_accession": "GCF_021172085.1_ASM2117208v1", 
  "length": 349, 
  "ccya": {
          "gr_accession": "NZ_CP089094.1", 
          "gr_start": 266112, 
          "gr_end": 267161, 
          "strand": "-", 
          "nucl_seq": "ATGGC....GCTAA", 
          "src": "Protein Homology", 
          "attr": "ID:cds-WP_084990073.1,Par...,product:cell envelope biogenesis protein OmpA,protein_id:WP_084990073.1,transl_table:11"
          }, 
   "calcyanin35": true, 
   "Observed iACC (y/n)": "y", 
   "calcyanin_flag": "Calcyanin with known N-ter", 
   "nter": "Z-type", 
   "cter": "Gly1,Gly2,Gly3", 
   "features": [
    {
      "feature_id": "nter", 
      "feature_type": "domain", 
      "start": 0, 
      "end": 113, 
      "source": 
      "blastp", 
      "score": 2.9e-79, 
      "ID": "Z-type", 
      "SRC": "NZ_CP020771.1_prot_WP_084990073.1_4271", 
      "coverage": 100
     },
     ], 
  "sequence": "MANPDSSAS....VS"}
```



[performance update : proceed by batch instead of assembly -> decrease the number of job, reduce queuing.]
# pyCALF
