#!/bin/bash

module load R/3.6.1-gimkl-2018b 
export R_LIBS=/scale_wlg_persistent/filesets/project/uoa02608/modules/GRIDSS_connector/R
export R_LIBS_USER=/scale_wlg_persistent/filesets/project/uoa02608/modules/GRIDSS_connector/R

#Rscript -e "install.packages('BiocManager', repos='http://cran.stat.auckland.ac.nz/') ; library('BiocManager') ; BiocManager::install(c('StructuralVariantAnnotation'))"
#Rscript -e "install.packages('/scale_wlg_persistent/filesets/project/uoa02608/modules/R_Lib_Arch/StructuralVariantAnnotation', repos = NULL, type = 'source')"
Rscript GRIDSS_connector.r $1 $2
