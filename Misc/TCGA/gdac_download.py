import requests
import time

cohorts = [ 'ACC','BLCA','BRCA','CESC','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','LAML','LGG','LIHC','LUAD','LUSC', 'OV', 'PAAD', 'PRAD','READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'UCEC', 'UCS' ]
destination_path = "D:/TCGA/"

for cohort in cohorts:

    file_set = {
##        'Clinical':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz'%(cohort, cohort),
##        'RawCounts':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
##        'NormCounts':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
##        'Mutations':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
##        'RPPA':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.RPPA_AnnotateWithGene.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
        'CNA1':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Merge_cna__cgh_1x1m_g4447a__mskcc_org__Level_3__segmentation_data_computation__seg.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
        'CNA2':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),
        'CNA3':'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz'%(cohort, cohort),}
    
    
    for k,f in file_set.items():
        response = requests.post( f, stream=True )
        if response.status_code==200:
            with open( destination_path + cohort + "/" + f.split("/")[-1], "wb" ) as q:
                for chunk in response.iter_content(chunk_size=512):
                    q.write( chunk )
        else:
            print( "Error: %s %s %d"%(cohort, k, response.status_code) )
        wait = time.sleep(1)
