#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#WGS dosage assessment script

#Read Input
bam=$1     #full path to indexed bam file
ID=$2      #string, ID for library
REF=$3     #string, specify either "h37" or "hg19"
check=$4   #string, either "summary", "genome", or "both"
OUTDIR=$5  #full path to output directory

# checks for appropriate input
if [ $# -eq 5 ]; then

#Hard-code variable paths
if [ ${REF} == "h37" ]; then
  int=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/dosage_intervals.bed
  dict=$DICT #/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.mainContigs.dict
else
  int=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/dosage_intervals.bed #/data/talkowski/Samples/SFARI/miscFiles/dosage_intervals.hg19.bed
  dict=$DICT #/data/talkowski/tools/ref/hg19_bwa07/hg19.lex.mainContigs.dict
fi

#Set params
h37=$REF #/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.fa
binsize=500000 #500kb bins by default
#TMPDIR=/scratch/miket/rlc47temp/tmp.files
if [ ${check} == "genome" ] || [ ${check} == "both" ]; then
  bins=`mktemp`
fi

#Create output directory
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
if ! [ -e ${OUTDIR}/${ID}_WGSdosageCheck ]; then
  mkdir ${OUTDIR}/${ID}_WGSdosageCheck
fi

#Get library size
total=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} )
sambamba view -f bam -F 'not secondary_alignment and not duplicate' ${bam} 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y> ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.bam
#different output for different specifications of "check"
if [ ${check} == "summary" ] || [ ${check} == "both" ]; then
  #Get counts
  intervals=`mktemp`
  cut -f1-4 ${int} > ${intervals}
  echo -e "STATUS [$(date)]: Getting counts"
  echo -e "NA\tNA\tNA\tTOTAL_LIBRARY\t${total}" > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed
  bedtools coverage -counts -sorted -b ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.bam -a ${intervals} >> ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed

  #Calculate fractions
  echo -e "STATUS [$(date)]: Calculating fractions"
  sed '1d' ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed | awk -v OFS="\t" -v total=${total} '{ print $1, $2, $3, $4, $5/total }' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.fractions.bed

  #Calculate summary fractions
  echo -e "STATUS [$(date)]: Calculating summary fractions"
  echo -e "IntervalClass\tFractionOfLibrary" > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.WGS_dosageSummary.txt
  for color in red blue yellow; do
    fgrep ${color} ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed | awk -v OFS="\t" -v total=${total} -v color=${color} '{ sum+=$5 } END { print color, sum/total }' >> ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.WGS_dosageSummary.txt
  done
fi
if [ ${check} == "genome" ] || [ ${check} == "both" ]; then
  #Bin genome
  contig_dir=`mktemp -d`
  while read contig length; do
    Rscript -e "options(scipen=1000); write.table(data.frame(as.character(\"${contig}\"),c(seq(0,${length},by=${binsize})),c(seq(${binsize},${length}+${binsize},by=${binsize}))),\"${contig_dir}/${contig}.bins\",row.names=F,col.names=F,sep=\"\t\",quote=F)"
    cat ${contig_dir}/${contig}.bins >> ${bins}
  done < <( awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk 'BEGIN{OFS="\t"};{ print $2, $4 }' )

  #Get counts
  echo -e "STATUS [$(date)]: Getting counts"
  exp_per_bin=$( echo "${total}/$( cat ${bins} | wc -l )" | bc )
  if [ ${REF} == "h37" ]; then
    bedtools coverage -sorted -counts -b ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.bam -a ${bins} | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed
  else
    bedtools coverage -sorted -counts -b ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.bam -a ${bins} | sed 's/^chr//g' | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' | sed 's/^/chr/g' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed
  fi

  #Calculate fractions
  echo -e "STATUS [$(date)]: Calculating fractions"
  awk -v OFS="\t" -v exp_per_bin=${exp_per_bin} '{ print $1, $2, $3, $4/exp_per_bin }' ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed

  #Clean up
  rm ${bins}
  rm -r ${contig_dir}
fi
rm ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.bam
#If input inappropriate, displays usage and exits
else
 echo -e "\n\nChromatin Dosage Anomaly Diagnostic Test\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
 echo "Usage:"
 echo "  WGScheckDosage.sh [path/to/bam] [ID] [h37/hg19] [summary/genome/both] [OUTDIR]"
 #echo ""
fi
