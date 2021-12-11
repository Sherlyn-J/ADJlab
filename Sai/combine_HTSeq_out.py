import os
import pandas as pd
import re

# combine htseq counts

def get_gene_list( fname="D:/gencode.v38.annotation.gtf" ):
    genes = []
    with open( fname ) as f:
        for line in f:
            if line.find( "HAVANA\tgene\t" ) != -1:
                d = line.split("\t")[-1].split( "; " )
                for item in d:
                    if item.find( "gene_name" )!= -1:
                        genes.append( item.split()[-1].replace('"','') )
	return genes


dirc = "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/HTSeq-count/"
fns  = os.listdir( "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/HTSeq-count/" )
genes= get_gene_list()
df   = pd.DataFrame( columns = [x.replace(".txt","") for x in fns] )

for fn in fns:
    col = fn.replace(".txt","")
    with open( dirc+fn, ) as q:
        for line in q:
            if line.find("processed.")==-1 and line.startswith("__")!=True:
                s = line.split()
                if s[0] in genes:
                    df.at[ s[0], col ] = s[1]


        
        

        

        
