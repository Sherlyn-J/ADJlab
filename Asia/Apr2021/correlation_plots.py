import pandas as pd
import numpy as np
import scipy.stats as sc
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="white")


# Load data
df = pd.read_csv( "D:/Michal-Asia/meta_annot.csv", header="infer", sep=",")
df.drop(columns="Sno", inplace=True)

# set subsets
o = ["D1","D2","D3","F1","F2","F3","F4","L1","L2","L3"]
P = ["MalignantB","MalignantB","MalignantB","MalignantB","MalignantB","MalignantB","MalignantB","HealthyB","HealthyB","HealthyB"]
origs = dict(zip( o, P ) )

# get spearman correlations
for k,v in origs.items():
    sub = df[ (df['orig.ident']==k) & (df.Population==v) ]
    corr = sc.spearmanr( sub.A3G, sub.CD20 )
    pl = sns.scatterplot( x="A3G", y="CD20", data=sub )
    pl.set_title( "%s: %s B-cells"%(sub.Sample.iat[0],v[:-1]) )
    # pl.legend(fontsize = 15, bbox_to_anchor= (1.03, 1), title="Delivery Type", title_fontsize = 18, shadow = True, facecolor = 'white')
    # plt.legend( [ "%d cells"%len(sub), "R=%0.2f\tp-val=%0.2e"%(corr[0], corr[1]) ] )
    plt.text(max(sub.A3G), max(sub.CD20)-(0.2*max(sub.CD20)), "%d cells\nSpearman R %0.2f\np-value %0.2e"%(len(sub),corr[0], corr[1]), horizontalalignment='right', size='medium')
    plt.savefig("D:/Michal-Asia/Scatter_A3G-CD20_%s.png"%sub.Sample.iat[0], dpi=360, bbox_inches='tight' )
    plt.clf()
    print( "%s: %s\tNCells: %d\tA3G:%0.2f-%0.2f\tCD20:%0.2f-%0.2f\tR=%0.2f\tp-val=%0.2e"%(k,v,len(sub),min(sub.A3G),max(sub.A3G),min(sub.CD20),max(sub.CD20), corr[0], corr[1] ) )

# DLBCL
sub = df[ (df['orig.ident'].isin(["D1","D2","D3"])) & (df.Population=="MalignantB") ]
corr = sc.spearmanr( sub.A3G, sub.CD20 )
pl = sns.scatterplot( x="A3G", y="CD20", data=sub )
pl.set_title( "DLBCL: Malignant B-cells" )
plt.text(max(sub.A3G), max(sub.CD20)-(0.2*max(sub.CD20)), "%d cells\nSpearman R %0.2f\np-value %0.2e"%(len(sub),corr[0], corr[1]), horizontalalignment='right', size='medium')
plt.savefig("D:/Michal-Asia/Scatter_A3G-CD20_DLBCL.png", dpi=360, bbox_inches='tight' )
plt.clf()

# FL
sub = df[ (df['orig.ident'].isin(["F1","F2","F3","F4"])) & (df.Population=="MalignantB") ]
corr = sc.spearmanr( sub.A3G, sub.CD20 )
pl = sns.scatterplot( x="A3G", y="CD20", data=sub )
pl.set_title( "FL: Malignant B-cells" )
plt.text(max(sub.A3G), max(sub.CD20)-(0.2*max(sub.CD20)), "%d cells\nSpearman R %0.2f\np-value %0.2e"%(len(sub),corr[0], corr[1]), horizontalalignment='right', size='medium')
plt.savefig("D:/Michal-Asia/Scatter_A3G-CD20_FL.png", dpi=360, bbox_inches='tight' )
plt.clf()

# rLN
sub = df[ (df['orig.ident'].isin(["L1","L2","L3"])) & (df.Population=="HealthyB") ]
corr = sc.spearmanr( sub.A3G, sub.CD20 )
pl = sns.scatterplot( x="A3G", y="CD20", data=sub )
pl.set_title( "rLN: Healthy B-cells" )
plt.text(max(sub.A3G), max(sub.CD20)-(0.2*max(sub.CD20)), "%d cells\nSpearman R %0.2f\np-value %0.2e"%(len(sub),corr[0], corr[1]), horizontalalignment='right', size='medium')
plt.savefig("D:/Michal-Asia/Scatter_A3G-CD20_rLN.png", dpi=360, bbox_inches='tight' )
plt.clf()
