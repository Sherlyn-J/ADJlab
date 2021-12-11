import pandas as pd

# get correlations, filter by abs(corr)>0.25 --> includes -0.25 and below as well
path = "C:/Users/csislyn/Dropbox (CSI (NUS))/Bryce/For Sherlyn/CCLE_RPPA_20181003.csv"
df = pd.read_csv( path )
x = df.corr()
x.to_csv( "C:/Users/csislyn/Dropbox (CSI (NUS))/Bryce/For Sherlyn/CCLE_RPPA_correlations.csv" )
y = x.loc[['Tuberin', 'Tuberin_pT1462', 'Akt_pT308', 'Akt_pS473', 'Akt'], ]
y = y.transpose()
y = y[ (abs(y.Akt)>0.25) | (abs(y.Akt_pS473)>0.25) | (abs(y.Akt_pT308)>0.25) | (abs(y.Tuberin)>0.25) | (abs(y.Tuberin_pT1462)>0.25) ]
y.to_csv( "C:/Users/csislyn/Dropbox (CSI (NUS))/Bryce/For Sherlyn/CCLE_RPPA_correlations_filtered.csv" )

# load mutations matrix into memory and map cell-line names to DepMap IDs
path_mutations = "C:/Users/csislyn/Dropbox (CSI (NUS))/Bryce/For Sherlyn/CCLE_mutations.csv"
path_mapfile = "C:/Users/csislyn/Dropbox (CSI (NUS))/Bryce/For Sherlyn/DepMap Line ID - For Mutation Mapping.csv"
mut_df = pd.read_csv( path_mutations )
map_df = pd.read_csv( path_mapfile )
out = mut_df.merge(map_df)
out.drop(columns="line",inplace=True)
