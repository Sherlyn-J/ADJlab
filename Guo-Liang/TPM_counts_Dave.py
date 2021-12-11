import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sc
import warnings
warnings.filterwarnings("ignore")

sns.set_context("paper")

# get gtf file and preprocess it for easy extraction of exon lengths
print( "Read in gtf" )
gtf = pd.read_csv( "C:/Users/csislyn/Downloads/hg38.ncbiRefSeq.gtf", sep="\t")
gtf = gtf[ gtf.type=='exon' ]
# insert column with gene name
print( "Add column for gene name" )
s   = []
for ndx, row in gtf.iterrows():
    s.append( str( row['desc'] ).split(";")[0].split()[1].replace('"','') )
gtf['gene'] = s
# insert column with gene length
print( "Add column for exon length" )
gtf['exon_length'] = abs(gtf.end-gtf.start)

raw_fn = 'C:/Users/csislyn/Dropbox (CSI (NUS))/Public-dataset-downloads/DLBCL/Dave-Cell-2017/Dave-Cell-2017_counts.csv'

# get raw counts and list of genes
print( "\tRead in main df" )
df  = pd.read_csv( raw_fn, header="infer", sep="," )
df.rename( columns={'gene_symbol':'Gene'}, inplace=True )
df  = df.drop_duplicates()
df  = df[~df.index.duplicated(keep='first')]
l   = list( df.Gene )
df  = df.set_index('Gene')

# calculate transcript length: sum of lengths of all exons
# main = { gn: sum(gtf[ gtf['gene']==gn ].exon_length) for gn in l }

# remove duplicates
print( "\tRemove duplicates" )
h = []
for ndx, row in df.iterrows():
    if row.name not in h:
        h.append( row.name )
    else:
        df.drop(index=row.name, inplace=True)

# get counts/transcript length for each gene
print( "\tdivide by transcript length" )
for ndx, row in df.iterrows():
    tcp = gtf[ gtf['gene']==row.name ]
    transcript_length = sum( tcp.exon_length ) if len(tcp) > 0 else 0
    if transcript_length > 0:
        df.loc[row.name] = pd.to_numeric(df.loc[row.name])/(transcript_length/1000)
    else:
        print( ndx, end=" " )
        df.loc[row.name] = pd.to_numeric(df.loc[row.name])/1.5 # average transcript length in kB

print( "\tscale by million" )
for col in df.columns:
    total = sum(df.loc[:,col])
    df.loc[:,col] = (df.loc[:,col]/total)*1000000
df.to_csv( "C:/Users/csislyn/Dropbox (CSI (NUS))/Public-dataset-downloads/DLBCL/Dave-Cell-2017/Dave-Cell-2017_TPMCounts.csv" )
