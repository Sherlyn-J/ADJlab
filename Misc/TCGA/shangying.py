# Process Shangying's data

import pandas as pd
import os
import numpy as np


print( "Load data..." )
file = "D:/TCGA/TCGA-26-Projects_combined_ASCAT_segmentation_purity-ploidy_0605-2021.txt"
df   = pd.read_csv( file, sep="\t", header="infer" )

print( "Calculate segment lengths..." )
df['Length'] = df.End - df.Start + 1

print( "Filter out segments of length 1..." )
df = df[ df.Length>1 ]
##df = df[ (df.project!="TCGA_OV") & (df.Chromosome!="chrY") ]

print( "Assigning CNA types..." )
df['X'] = df.apply( lambda row: row.barcode[:12], axis=1 )
##df['CNAtype']= df.apply( lambda row: 'GAIN' if (np.log2(row.Copy_Number/row.Ploidy) > 1) else 'LOSS' if (np.log2(row.Copy_Number/row.Ploidy) < -1) else 'NEUTRAL', axis = 1 )
df['CNAtype']= df.apply( lambda row: 'GAIN' if (row.Copy_Number-row.Ploidy>0.6) else 'LOSS' if (row.Copy_Number-row.Ploidy<-0.6) else 'NEUTRAL', axis = 1 )

print( "Assigning LOH..." )
df['LOH'] = df.apply( lambda row: 'TRUE' if row.Minor_Copy_Number==0 else "FALSE", axis=1 )

##for ndx, row in df.iterrows():
##    if row.Copy_Number - row.Ploidy > 0.6:
##        row.CNAtype = 'GAIN'
##    elif row.Copy_Number - row.Ploidy < -0.6:
##        row.CNAtype = 'LOSS'

cancerTypes = [ x for x in os.listdir(".") if os.path.isdir( os.path.join( ".",x ) ) and x!='Aggregate' ]

print( "Partition data..." )
for cancer in cancerTypes:
    sub = df[ df.project.str.contains( cancer ) ]
    sub.to_csv( "D:/TCGA/%s/%s_shangying.csv"%(cancer,cancer), sep="\t", index=False )

print( "Done" )
