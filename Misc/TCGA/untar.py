import gzip
import zipfile
import tarfile
import os
import shutil

working_path = "D:/TCGA/"
cohorts = [ 'ACC','BLCA','BRCA','CESC','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','LAML','LGG','LIHC','LUAD','LUSC', 'OV', 'PAAD', 'PRAD','READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'UCEC', 'UCS' ]

# unzip gz files to tar
##for cohort in cohorts:    
##    directory = working_path + cohort + "/"
##    for file in os.listdir( path=directory ):
##        with gzip.open(directory + file, 'rb') as f_in:
##            with open(directory + file[:-3], 'wb') as f_out:
##                shutil.copyfileobj(f_in, f_out)
##
# uncompress tar files
for cohort in cohorts:
    directory = working_path + cohort + "/"
    for file in [fn for fn in os.listdir( path=directory ) if fn.endswith(".tar")]:
        if tarfile.is_tarfile(directory+file):
            with tarfile.open(directory+file, 'r') as f_in:
                for m in f_in.getmembers():
                    path = directory + m.path.split("/")[-1]
                    if m.isdir():
                        tarfile.TarFile.makedir( f_in, directory + m.path.split("/")[-1] )
                        path = directory + m.path.split("/")[-1]
                    else:
                        try:
                            tarfile.TarFile.extract( f_in, m, path=directory )
                        except:
                            tarfile.TarFile.extract( f_in, m )
