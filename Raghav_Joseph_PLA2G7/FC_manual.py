import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sc
from bioinfokit import analys, visuz

sns.set_theme(style="whitegrid")

# import methods and data
# get r51 quartile assignments
# filter out immune genes
# filter out q1 and Q4 r51
# calculate median fold change: log2( medianq4/medianq1 )
# calculate wilcoxon p-value
# apply bonferroni correction
# set up data for volcano plot in R

def assign_quartile( array ):
    cuts = np.quantile(r51, [ .25, .5, .75,])
    assigned = []
    for i in array:
        if i > cuts[-1]:
            assigned.append(4)
        elif i > cuts[-2]:
            assigned.append(3)
        elif i > cuts[-3]:
            assigned.append(2)
        else:
            assigned.append(1)
    return assigned


df = pd.read_csv( "D:/Michal_Jan2021/SAMIT/samitnano.csv" )
immune = list(df[ df.Annotation=="Immune" ]["Probe Name"])
df_ = df.set_index("Probe Name").transpose()

df = df[ df.Annotation=="Immune" ]
df.drop('Annotation', axis=1, inplace=True)

df_.drop('Annotation', axis=0, inplace=True)
r51 = list(df_.RAD51)
df_['r51q'] = assign_quartile( list(df_.RAD51) )


low_samples  = df_[ df_.r51q==1 ].index
high_samples = df_[ df_.r51q==4 ].index

results = pd.DataFrame( columns=["Probe Name", "MedianQ4", "MedianQ1", "log2FC", "pval",] )

for ndx, gene in enumerate(immune):
    results.loc[ndx] = [ gene,
                      df_[ df_.r51q==4 ][gene].median(),
                      df_[ df_.r51q==1 ][gene].median(),
                      np.log2(df_[ df_.r51q==4 ][gene].median()/df_[ df_.r51q==1 ][gene].median()),
                      sc.wilcoxon( list(df_[ df_.r51q==4 ][gene]), list(df_[ df_.r51q==1 ][gene]) ).pvalue ]

print( results )
signames = results[ (results.pval < 0.0006) & ((results.log2FC > 1)|(results.log2FC <1)) ]['Probe Name']

#visuz.gene_exp.volcano(df=results, lfc='log2FC', pv='pval', color=("#077E97", "grey", "#94641F"),
#                       genenames=list(signames) )

results['log10(pvalue)'] = -1*np.log10(results.pval)
results['select']        = results.apply( lambda row: "Q4 high" if row.log2FC>.5 and row.pval<.0006 else "Q4 low" if row.log2FC<-.5 and row.pval<.0006 else "", axis=1)

ax = sns.scatterplot(data=results, x="log2FC", y="log10(pvalue)", hue="select", palette={"Q4 low":"#077E97","":"grey","Q4 high":"#94641F"}, alpha=0.7)
ax.grid(b=None)
plt.axvline(x=0.5, color="grey", linestyle="--",lw=0.5)
plt.axvline(x=-0.5, color="grey", linestyle="--",lw=0.5)
plt.axhline(y=-1*np.log10(0.0006), color="grey", linestyle="--",lw=0.5)
plt.ylim(0,20)
plt.legend(loc=7)
pdx = (ax.get_xlim()[1]-ax.get_xlim()[0])/100
pdy = (ax.get_ylim()[1]-ax.get_ylim()[0])/50
for ndx, row in results.iterrows():
    if row.select=="Q4 high":
        ax.text(row.log2FC+(pdx*0.6), row['log10(pvalue)']+(pdy*0.1), row["Probe Name"], color="#94641F", fontdict={'size':7.5})
    elif row.select=="Q4 low":
        ax.text(row.log2FC+(pdx*0.6), row['log10(pvalue)'], row["Probe Name"], color="#077E97", fontdict={'size':7.5})
plt.savefig("D:/Michal_Jan2021/SAMIT/samit_immune_volcano.png", dpi=400)
