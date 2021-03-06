#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

# Script to create bedfile of all exons from transcripts which overlap input region

##NOTE: REQUIRES MANUAL INSTALLATION OF THE R LIBRARY "PLYR"

##INPUT##
ref=$1          ##ref flat; either hg18 or hg19
chr=$2          ##chromosome
start=$3        ##start position
stop=$4         ##stop position
OUTFILE=$5      ##output bed file

# checks for appropriate input
if [ $# -eq 5 ]; then

#Loads appropriate version of R
#module load R/R-3.0.0

##Creates tmp file
TMPFILE=`mktemp`
#echo "${TMPFILE}"

##Writes all overlapping entries from RefFlat to tmp file
if [ ${ref} == "hg18" ]; then
 awk -v chr=${chr} -v A=${start} -v B=${stop} \
 'BEGIN{FS=","}; $3 ~ /'\"$chr\"'/ { if ($5<=B && $6>=A) print}' \
 /data/talkowski/TGDB/BACKUP/hg18_refFlat_3_13_2014.csv > ${TMPFILE}
elif [ ${ref} == "hg19" ]; then
 awk -v chr=${chr} -v A=${start} -v B=${stop} -v FS="\t" '{ if ($2<=B && $3>=A && $1==$chr) print}' ${refFlat} > ${TMPFILE}
else
 echo "ERROR: PLEASE SPECIFY EITHER hg18 OR hg19"
fi

##Executes Rscript to build & write exon table
Rscript -e "library(plyr); if(file.info(\"${TMPFILE}\")\$size > 0){tmp <- read.table(\"${TMPFILE}\",header=F,stringsAsFactors=F); \
tmp_table = unique(rbind.fill(suppressWarnings(apply(tmp,1,function(vals){unique(data.frame(vals[1],unlist(strsplit(vals[11],split=\",\")),unlist(strsplit(vals[12],split=\",\")),vals[4]))})))); \
write.table(tmp_table,\"${OUTFILE}\",sep=\"\t\",quote=F,row.names=F,col.names=F)}"

##Cleans tmpfile
rm -rf ${TMPFILE}

# if inappropriate input; displays useage
else
 echo -e "\n\nExon Builder\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
 echo "Usage:"
 echo "  build_exons.sh [hg18/hg19] [chr] [start] [stop] [OUTFILE]"
 echo ""
fi
