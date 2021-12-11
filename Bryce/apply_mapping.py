import pandas as pd

df = pd.read_csv( "D:/Bryce/Everolimus-sensitivity/GSE155923_SCLC_BETi_mTORi_combo_Chen2020_normed_counts.tsv",sep="\t",header="infer")

ensembl_to_hgnc = {}
with open( "D:/Bryce/Everolimus-sensitivity/HGNC_GS.txt", ) as f:
    for line in f:
        s = line.rstrip().split("\t")
        try:
            if s[1] != '':
                ensembl_to_hgnc[ s[0] ] = s[1]
        except:
            a = next

df['HGNC'] = df.apply( lambda row: ensembl_to_hgnc[ row['ensgene'] ] if row['ensgene'] in ensembl_to_hgnc else '' , axis=1 )
df.to_csv( "D:/Bryce/Everolimus-sensitivity/GSE155923_HGNC_mapped.txt", sep="\t" )
