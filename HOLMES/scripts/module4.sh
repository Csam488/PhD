#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Module 4 (Physical Depth CNV Calling)

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}
cd ${WRKDIR}

#module rm gcc/4.9.0
#module rm gcc-4.4
#module load gcc/4.9.0

#Load cnMOPS modules
#module unload R/3.1.0
#module load R/R-3.0.0

#Submit cnMOPS - autosomes
echo -e "STATUS [$(date)]: Submitting cnMOPS for autosomes..."
#Rscript -e "source('https://bioconductor.org/biocLite.R') ; biocLite('cn.mops') ; library(cn.mops)"
#Rscript -e "install.packages('BiocManager', repos='http://cran.stat.auckland.ac.nz/') ; library('BiocManager') ; chooseBioCmirror() ; BiocManager::install(c('cn.mops'))"
#Install cn.mops
#Rscript -e "install.packages('${liWGS_SV}/R_Lib_Arch/zlibbioc', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/BiocGenerics', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/BiocParallel', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/Biobase', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/S4Vectors', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/IRanges', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/XVector', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/Rhtslib', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/Biostrings', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/GenomeInfoDbData', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/GenomeInfoDb', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/GenomicRanges', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/Rsamtools', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/exomeCopy', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/cn.MOPS', repos = NULL, type = 'source')"
#Install rtracklayer
#Rscript -e "install.packages('${liWGS_SV}/R_Lib_Arch/DelayedArray', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/SummarizedExperiment', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/GenomicAlignments', repos = NULL, type = 'source');install.packages('${liWGS_SV}/R_Lib_Arch/rtracklayer', repos = NULL, type = 'source')"
mkdir ${WRKDIR}/cnMOPS
cov=${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed
for contig in $( seq 1 22 ); do #<=============================================
  mkdir ${WRKDIR}/cnMOPS/${contig}
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed
  for binsize in 1 3 10 30; do
    echo -e "STATUS [$(date)]: For bin size ${binsize}kb..."
    Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID} ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed ${WRKDIR}/cnMOPS/${contig}/
  done
done

#Submit cnMOPS - allosomes
echo -e "STATUS [$(date)]: Submitting cnMOPS for allosomes..."
if [ ${other_assign} == "MALE" ]; then
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="M" || $14=="O") print $1 }' > ${WRKDIR}/males.list
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="F") print $1 }' > ${WRKDIR}/females.list
else
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="M") print $1 }' > ${WRKDIR}/males.list
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="F" || $14=="O") print $1 }' > ${WRKDIR}/females.list
fi
echo $PWD
Midx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/males.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )
Fidx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/females.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )
for contig in X Y; do
  rm -rf ${WRKDIR}/cnMOPS/${contig}
  mkdir ${WRKDIR}/cnMOPS/${contig}
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) | cut -f ${Midx} > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.M.bed
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) | cut -f ${Fidx} > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.F.bed
  for binsize in 1 3 10 30; do
    echo -e "STATUS [$(date)]: For bin size ${binsize}kb..."
    Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID}_M ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.M.bed ${WRKDIR}/cnMOPS/${contig}/
    Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID}_F ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.F.bed ${WRKDIR}/cnMOPS/${contig}/
  done
done

#Gate until complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cnMOPS" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cnMOPS" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo "$(date): INCOMPLETE"
#    echo "$(date): Waiting on ${GATEcount} jobs to complete"
#    GATEwait=0
#  fi
#done

#Merge cnMOPS calls
echo -e "STATUS [$(date)]: Merging..."
for contig in $( seq 1 22 ); do
  for binsize in 1 3 10 30; do
    sed -i 's/;/\t/g' ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.${binsize}000bpBins.cnMOPS.gff
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
  done
done
for contig in X Y; do
  for binsize in 1 3 10 30; do
    sed -i 's/;/\t/g' ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_M.${binsize}000bpBins.cnMOPS.gff
    sed -i 's/;/\t/g' ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_F.${binsize}000bpBins.cnMOPS.gff
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_M.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_F.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
  done
done
#Split cnMOPS Calls by sample
echo -e "STATUS [$(date)]: Splitting by sample..."
mkdir ${WRKDIR}/cnMOPS/cnMOPS_calls
while read bam ID; do
  mkdir ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}
  fgrep -w "${ID}" ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff | grep 'CN[0-1]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed
  fgrep -w "${ID}" ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff | grep 'CN[3-9]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed
done < ${WRKDIR}/cnMOPS.input

echo -e "$(date)" > ${OUTDIR}/checkpoints/module4.done