import pandas as pd
import warnings
warnings.filterwarnings("ignore")

working_path = "D:/TCGA/" # 'ACC','DLBC','LAML','UCS', 
cohorts = [ 'BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC','KICH','KIRC','LGG','LIHC','LUAD','LUSC', 'OV', 'PAAD', 'PRAD','READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'UCEC' ]

for cohort in cohorts:
    print( cohort, end="... " )
    clinical_fn = 'D:\TCGA\%s\gdac.broadinstitute.org_%s.Clinical_Pick_Tier1.Level_4.2016012800.0.0\All_CDEs.txt'%( cohort, cohort )
    raw_fn = 'D:\TCGA\%s\gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0\%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt'%( cohort, cohort, cohort )

    df = pd.read_table(raw_fn, header="infer")
    cd = pd.read_csv( clinical_fn, header="infer", sep="\t" )

    ## Format main data
    df.set_index( "Hybridization REF", inplace=True )
    df = df.transpose()
    df = df[ df.gene_id=="raw_count" ]
    df.index = [ "-".join( x.split("-")[:3] ) for x in df.index ]
    df.drop( columns='gene_id', inplace=True )
    df.columns = [ x.split("|")[0] for x in df.columns ]
    df.drop( columns='?', inplace=True )
    df = df.transpose()

    ## standardize patient names
    df.columns = [ x[:12] for x in df.columns ]
    df = df.transpose()
    df = df[~df.index.duplicated(keep='first')]
    df = df.transpose()

    ## set index
    cd.set_index( 'bcr_patient_barcode', inplace=True )
    cd.columns = [ x.upper() for x in cd.columns ]
    
    ## Filter by drop, reorder to match clindata
    patient_list = [x.upper() for x in list(cd.columns)] # or patient_id
    newdf = pd.DataFrame( df, columns=[x for x in df.columns if x in cd.columns] )
    newcd = pd.DataFrame( cd, columns=[x for x in cd.columns if x in df.columns] )        

##    ## rearrange according to clinical data
##    newdf = pd.DataFrame( df, columns=cd.columns )
##    newdf.dropna( axis=1, inplace=True )

##    ## filter clinical data
##    newcd = pd.DataFrame( cd, columns=newdf.columns[1:] )

    # write
    newdf.to_csv( working_path + "/%s/%s_RawCounts.csv"%(cohort,cohort),sep="\t" )
    newcd.to_csv( working_path + "/%s/%s_Clinical.csv"%(cohort,cohort),sep="\t" )
    print( "Done" )
