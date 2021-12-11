import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from adjustText import adjust_text

sns.set(style = 'whitegrid')

pathways_to_mark = [ 'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION' ] # 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2', 'HALLMARK_TNFA_SIGNALLING_VIA_NFKB',

def preformat( file ):
    with open( file ) as q:
        data = q.readlines()
    data = [ line.replace( "\tDetails ...\t", "\t" ) for line in data ]
    data = [ line.replace( "\t\t", "\t" ) for line in data ]
    data = [ line.replace( "\tGS<br> follow link to MSigDB\t", "\t" ) for line in data ]
    with open( file, 'w' ) as q:
        q.write( "".join(data) )

tname = "OVCAR3_KO_untreated"
cname = "OVCAR3_P_untreated"
        
ttmt_file = "D:\Sherlyn-bioinfo\mitochondria\RNA-seq_processed\DEGs_29oct21\OVCAR3_KO_vs_P.Gsea.1635478252344\gsea_report_for_OVCAR3_KO_untreated_1635478252344.tsv"
ctrl_file = "D:\Sherlyn-bioinfo\mitochondria\RNA-seq_processed\DEGs_29oct21\OVCAR3_KO_vs_P.Gsea.1635478252344\gsea_report_for_OVCAR3_P_untreated_1635478252344.tsv"

preformat( ttmt_file )
preformat( ctrl_file )
ttmt = pd.read_csv( ttmt_file, sep="\t", header="infer" )
ctrl = pd.read_csv( ctrl_file, sep="\t", header="infer" )

x1 = list(ttmt.NES)
x2 = list(ctrl.NES)
y  = list( -1*np.log10(ttmt['FDR q-val']) ) + list(-1*np.log10(ctrl['FDR q-val']))

# plt.scatter(x = ttmt.NES, y = -1*np.log10(ttmt['FDR q-val']), color = 'grey', alpha = 0.4)
# plt.scatter(x = ctrl.NES, y = -1*np.log10(ctrl['FDR q-val']), color = 'grey', alpha = 0.4)
cmap = mpl.cm.viridis
norm = plt.Normalize(vmin=min(x2),vmax=max(x1))
plt.scatter(x = x1+x2, y = y, c=x1+x2, norm=norm, cmap=cmap, alpha=0.4)

texts = []

# Label selected points: iterate both dfs
for ndx, row in ttmt.iterrows():
    if row.NAME in pathways_to_mark:
        texts.append( plt.text( row.NES, -1*np.log10(row['FDR q-val']), row.NAME.replace("HALLMARK_",""), fontsize=7 ) )

for ndx, row in ctrl.iterrows():
    if row.NAME in pathways_to_mark:
        texts.append( plt.text( row.NES, -1*np.log10(row['FDR q-val']), row.NAME.replace("HALLMARK_",""), fontsize=7 ) )

# title and axes labels
plt.title( "Hallmark pathway enrichment: %s vs %s"%(tname, cname) )
plt.xlabel( "Normalized ES" )
plt.ylabel( "-log10 (FDR q-value)" )
plt.annotate(tname, ((max(x1)/2),max(y)-0.05),textcoords="offset points", xytext=(0,-2), ha='center', fontweight='bold')
plt.annotate(cname, ((min(x2)/2)-0.2,max(y)-0.05),textcoords="offset points", xytext=(0,-2), ha='center', fontweight='bold')
plt.xlim( (-1.5,1.5) )
# vertical line
plt.axvline(color="grey")
# set grid None
plt.grid(None)
adjust_text(texts, only_move={'points':'', 'texts':'x,y'}, arrowprops=dict(arrowstyle="-", color='r', lw=0.5))
plt.savefig( "D:\Sherlyn-bioinfo\mitochondria\RNA-seq_processed\DEGs_29oct21\%s_%s-hm-pathways_2.png"%(tname,cname), dpi=400 )
