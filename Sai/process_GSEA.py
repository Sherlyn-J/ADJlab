import pandas as pd
import numpy as np

wd = "D:/Michal_Jul2021/sc_analysis/GSEA/"
names = {}
with open( wd+"hallmarks.txt", ) as f:
    for line in f:
        s = line.split("\t")
        names[ s[0] ] = s[1].rstrip()

print( names )

def preformat( file ):
    with open( file ) as q:
        data = q.readlines()
    data = [ line.replace( "\tDetails ...\t", "\t" ) for line in data ]
    data = [ line.replace( "\t\t", "\t" ) for line in data ]
    data = [ line.replace( "\tGS<br> follow link to MSigDB\t", "\t" ) for line in data ]
    with open( file, 'w' ) as q:
        q.write( "".join(data) )


# record ttmt, cntrl, ttmt_file, ctrl_file for each comparison made
experiments = {'DLBC M+2+6- vs all':{'tname':'M+2+6-',
                                  'cname':'all',
                                  'ttmt_file':wd+'DLBC_M26_vs_all.GseaPreranked.1632291305794/gsea_report_for_na_pos_1632291305794.tsv',
                                  'ctrl_file':wd+'DLBC_M26_vs_all.GseaPreranked.1632291305794/gsea_report_for_na_neg_1632291305794.tsv',},
               'DLBC M+2+6- vs M+':{'tname':'M+2+6-',
                                  'cname':'M+',
                                  'ttmt_file':wd+'DLBC_M26_vs_M.GseaPreranked.1632291342370/gsea_report_for_na_pos_1632291342370.tsv',
                                  'ctrl_file':wd+'DLBC_M26_vs_M.GseaPreranked.1632291342370/gsea_report_for_na_neg_1632291342370.tsv',},
               'DLBC M+2+6- vs 2+':{'tname':'M+2+6-',
                                  'cname':'2+',
                                  'ttmt_file':wd+'DLBC_M26_vs_2.GseaPreranked.1632291353514/gsea_report_for_na_pos_1632291353514.tsv',
                                  'ctrl_file':wd+'DLBC_M26_vs_2.GseaPreranked.1632291353514/gsea_report_for_na_neg_1632291353514.tsv',},
                'rLN M+2+6- vs all':{'tname':'M+2+6-',
                                  'cname':'all',
                                  'ttmt_file':wd+'rLN_M26_vs_all.GseaPreranked.1632291372606/gsea_report_for_na_pos_1632291372606.tsv',
                                  'ctrl_file':wd+'rLN_M26_vs_all.GseaPreranked.1632291372606/gsea_report_for_na_neg_1632291372606.tsv',},
                'rLN M+2+6- vs M+':{'tname':'M+2+6-',
                                  'cname':'M+',
                                  'ttmt_file':wd+'rLN_M26_vs_M.GseaPreranked.1632291388561/gsea_report_for_na_pos_1632291388561.tsv',
                                  'ctrl_file':wd+'rLN_M26_vs_M.GseaPreranked.1632291388561/gsea_report_for_na_neg_1632291388561.tsv',},
                'rLN M+2+6- vs 2+':{'tname':'M+2+6-',
                                  'cname':'2+',
                                  'ttmt_file':wd+'rLN_M26_vs_2.GseaPreranked.1632291398196/gsea_report_for_na_pos_1632291398196.tsv',
                                  'ctrl_file':wd+'rLN_M26_vs_2.GseaPreranked.1632291398196/gsea_report_for_na_neg_1632291398196.tsv',}, }

for en, ed in experiments.items():

    preformat( ed['ttmt_file'] )
    preformat( ed['ctrl_file'] )
    ttmt = pd.read_csv( ed['ttmt_file'], sep="\t", header="infer" )
    ctrl = pd.read_csv( ed['ctrl_file'], sep="\t", header="infer" )
    new = pd.concat( [ttmt, ctrl], axis=0 )
    new['NAME2'] = new.apply( lambda row: names[row.NAME], axis=1 )
    new.to_csv( wd+en+"_gsea_processed.csv", sep="\t", index=False )

    

