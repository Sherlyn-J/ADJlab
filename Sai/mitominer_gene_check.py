import pandas as pd
from Bio import SeqIO


##df = pd.read_csv( "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/3433_merged_isoform_TPM.txt", sep="\t", header="infer" )
##df['Isoform2'] = df.apply( lambda row: row.Isoform.split("|")[0], axis=1 )
##genes_to_search = [ 'TRNF', 'MT-TF', 'RNR1', 'MTRNR1', 'TRNV', 'MTTV', 'RNR2', 'MTRNR2', 'TRNL1', 'MTTL1', 'ND1', 'MTND1', 'TRNI', 'MTTI', 'TRNQ', 'MTTQ', 'TRNM', 'MTTM', 'ND2', 'MTND2', 'TRNW', 'MTTW', 'TRNA', 'MTTA', 'TRNN', 'MTTN', 'TRNC', 'MTTC', 'TRNY', 'MTTY', 'COX1', 'COI', 'MTCO1', 'TRNS1', 'MTTS1', 'TRND', 'MTTD', 'COX2', 'COII', 'MTCO2', 'TRNK', 'MTTK', 'ATP8', 'MTATP8', 'ATPASE8', 'ATP6', 'ATPASE6', 'MTATP6', 'COX3', 'COIII', 'MTCO3', 'TRNG', 'MTTG', 'ND3', 'MTND3', 'TRNR', 'MTTR', 'ND4L', 'MTND4L', 'ND4', 'MTND4', 'TRNH', 'MTTH', 'TRNS2', 'MTTS2', 'TRNL2', 'MTTL2', 'ND5', 'MTND5', 'ND6', 'MTND6', 'TRNE', 'MTTE', 'CYTB', 'MTCYB', 'TRNT', 'MTTT', 'TRNP', 'MTTP', ]
##
##print( df[ df.Gene.isin( genes_to_search ) ] )
##
##new_genes_ids = []
##for record in SeqIO.parse( "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/NC012920.1.gb", "gb"):
##    for feature in record.features:
##        if feature.type=="gene":
##            new_genes_ids += feature.qualifiers['gene']
##            if 'gene_synonym' in feature.qualifiers:
##                new_genes_ids += feature.qualifiers['gene_synonym']


##mitominer = pd.read_csv( "D:/Sherlyn-bioinfo/mitochondria/mitominer.tsv", sep="\t", header="infer" )
##mitominer_MT = mitominer[ mitominer.Chromosome=="MT" ]
##genes = list( mitominer['Gene Symbol'] )
##print( df[ df.Isoform2.isin( genes ) ] )


# counts, TPMs
df_counts = pd.read_csv( "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/3433_merged_gene_counts.txt", sep="\t", header="infer" )
df_tpms   = pd.read_csv( "D:/Sherlyn-bioinfo/mitochondria/RNA-seq_processed/3433_merged_isoform_TPM.txt", sep="\t", header="infer" )
df_tpms['Isoform2'] = df_tpms.apply( lambda row: row.Isoform.split("|")[0], axis=1 )

# mitominer
mitominer = pd.read_csv( "D:/Sherlyn-bioinfo/mitochondria/mitominer.tsv", sep="\t", header="infer" )

# 
merged = pd.merge( mitominer, df_tpms, left_on="Gene Symbol", right_on="Isoform2" )
ids = merged[ ['Ensembl Primary Identifier', 'Gene Symbol', 'Chromosome', 'NCBI Gene ID'] ]
ids.drop_duplicates(inplace=True)

df_totals = pd.merge( ids, merged.groupby( "Ensembl Primary Identifier" ).sum(), on="Ensembl Primary Identifier" )
df_totals[ ['Gene Symbol', 'Ges1_A_10_E1', 'Ges1_A_10_E2', 'Ges1_A_D_E1', 'Ges1_A_D_E2',
            'GES1_P_10_E1', 'GES1_P_10_E2', 'GES1_P_D_E1', 'GES1_P_D_E2',
            'MCF10A_AIAKO_2', 'MCF10A_KO_1', 'MCF10A_Parent_2',
            'OVCAR3_A1AKO_1', 'OVCAR3_A1AKO_2', 'OVCAR3_P_1', 'OVCAR3_P_2'] ].to_csv( "D:/Sherlyn-bioinfo/mitochondria/MitoMiner_results.txt", sep="\t", index=False )
