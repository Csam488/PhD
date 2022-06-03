#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Calculate either nucleotide or physical coverage for a list of input bams

#Reads input
list=$1          #list of bams for analysis.  Note that this should be tab-delimmed 
                 #with two columns: first col = full path, second col = sample ID
dict=$2          #ref dict
mode=$3          #either "nucleotide" or "physical"
binsize=$4       #bin size for genome segmentation
ID=$5            #group ID
OUTDIR=$6		     #output directory

# checks for appropriate input
if [ $# -eq 6 ]; then

  #loads modules
  #module load BEDTools_2.17
  #module load samtools-0.1.14

  #creates directory tree
  WORKING=`mktemp -d -p ${TMPDIR}`
  echo -e "\nWriting temporary files to ${WORKING}"
  if [ ! -d "$OUTDIR" ]; then
    mkdir ${OUTDIR}
  fi
  if [ ! -d "${OUTDIR}/raw_coverages" ]; then
    mkdir ${OUTDIR}/raw_coverages
  fi
  if [ ! -d "${OUTDIR}/logfiles" ]; then
    mkdir ${OUTDIR}/logfiles
  fi
  echo -e "Writing output files to ${OUTDIR}"

  #creates contig list
  awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk \
  'BEGIN{OFS="\t"};{ print $2, $4 }' > ${WORKING}/contig.list
  num_contigs=$( cat ${WORKING}/contig.list | wc -l )

  #splits proper primary alignment pairs from each input bam
  if [ ${mode} == "physical" ]; then
    echo -e "EXTRACTING PROPER PAIRS..."
    while read bam sample; do
      mkdir ${WORKING}/${sample}
      echo "${sample}"
      while read contig_line; do
        contig=$( echo ${contig_line} | awk '{ print $1}' )
        sambamba view -f bam -F \
        'paired and proper_pair and not (unmapped or mate_is_unmapped or duplicate or secondary_alignment)' \
        ${bam} ${contig} | samtools sort \
        -T ${WORKING}/${sample}.${contig}.screened -O bam /dev/stdin > \
        ${WORKING}/${sample}/${sample}.${contig}.screened.bam
      done < ${WORKING}/contig.list
    done < ${list}
  fi

  #get coverage for nucleotide mode
  if [ ${mode} == "nucleotide" ]; then
    echo -e "\nCALCULATING COVERAGE PER CONTIG...\n"
    while read contig length; do
      while read bam sample; do
        sambamba view -f bam -F 'not secondary_alignment and not duplicate' ${bam} ${contig} | \
        bedtools coverage -counts -b - -a ${WORKING}/${contig}.bins | sort -nk2,2 > \
        ${WORKING}/${sample}.${contig}.coverage.bed
      done < ${list}
    done < ${WORKING}/contig.list
  fi

  #creates segmentation bed
  echo -e "\nSEGMENTING REFERENCE...\n"
  awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk 'BEGIN{OFS="\t"};{ print $2, $4 }' > ${WORKING}/contig.list
  while read contig length; do
   echo "BINNING ${contig}"
   Rscript -e "options(scipen=1000); write.table(data.frame(as.character(\"${contig}\"),c(seq(0,${length},by=${binsize})),c(seq(${binsize},${length}+${binsize},by=${binsize}))),\"${WORKING}/${contig}.bins\",row.names=F,col.names=F,sep=\"\t\",quote=F)"
   cat ${WORKING}/${contig}.bins >> ${WORKING}/binned_genome.bed
  done < ${WORKING}/contig.list

  ##Gate (20 second check; 5 minute report)
  #echo -e "\n"
  #GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_step1" | wc -l)
  #GATEwait=0
  #until [[ $GATEcount == 0 ]]; do
  # sleep 20s
  # GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_step1" | wc -l)
  # GATEwait=$[$GATEwait +1]
  # if [[ $GATEwait == 15 ]]; then
  #  echo "$(date): INCOMPLETE"
  #  echo "$(date): Waiting on ${GATEcount} step 1 jobs to complete"
  #  GATEwait=0
  # fi
  #done

  #merges filtered alignments split by contig (if necessary), otherwise sorts and renames bams
  if [ ${mode} == "physical" ]; then
    if [ $num_contigs -gt 1 ]; then 
     echo -e "\nMERGING & SORTING FILTERED ALIGNMENTS...\n"
     while read bam sample; do
      sambamba merge -l 6 ${WORKING}/${sample}/${sample}.merged.bam ${WORKING}/${sample}/${sample}.*.screened.bam
      sambamba sort -n --tmpdir=${WORKING} -l 6 -o ${WORKING}/${sample}.screened.nsort.bam ${WORKING}/${sample}/${sample}.merged.bam
     done < ${list}
    elif [ $num_contigs -eq 1 ]; then
     echo -e "\nSORTING FILTERED ALIGNMENTS...\n"
     contig=$( awk '{ print $1 }' ${WORKING}/contig.list )
     while read bam sample; do
      sambamba sort -m 4GB -n -l 6 -o ${WORKING}/${sample}.screened.nsort.bam ${WORKING}/${sample}/${sample}.${contig}.screened.bam
     done < ${list}
    fi
    ##Gate (20 second check; 5 minute report)
    #echo -e "\n"
    #GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_mergefiltered" | wc -l)
    #GATEwait=0
    #until [[ $GATEcount == 0 ]]; do
    # sleep 20s
    # GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_mergefiltered" | wc -l)
    # GATEwait=$[$GATEwait +1]
    # if [[ $GATEwait == 15 ]]; then
    #  echo "$(date): INCOMPLETE"
    #  echo "$(date): Waiting on ${GATEcount} alignment mergers to complete"
    #  GATEwait=0
    # fi
    #done
  fi

  #creates fragments
  if [ ${mode} == "physical" ]; then
    echo -e "\nSYNTHESIZING FRAGMENTS...\n"
    while read bam sample; do
      bedtools bamtobed -bedpe -i ${WORKING}/${sample}.screened.nsort.bam | \
      awk -v OFS="\t" '{ if ($2>=0 && $6>$2) print $1, $2, $6 }' > \
      ${WORKING}/${sample}.fragments.bed
    done < ${list}
    ##Gate (20 second check; 5 minute report)
    #echo -e "\n"
    #GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_synthesizefragments" | wc -l)
    #GATEwait=0
    #until [[ $GATEcount == 0 ]]; do
    # sleep 20s
    # GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_synthesizefragments" | wc -l)
    # GATEwait=$[$GATEwait +1]
    # if [[ $GATEwait == 15 ]]; then
    #  echo "$(date): INCOMPLETE"
    #  echo "$(date): Waiting on ${GATEcount} fragment synthesis jobs to complete"
    #  GATEwait=0
    # fi
    #done
  fi

  #calculates coverage
  if [ ${mode} == "physical" ]; then
    echo -e "\nCALCULATING COVERAGE...\n"
    while read bam sample; do
     bedtools coverage -counts -b ${WORKING}/${sample}.fragments.bed \
     -a ${WORKING}/binned_genome.bed > ${WORKING}/${sample}.unsorted.coverage.bed; \
     sort -V -k1,1 -k2,2n -k3,3n ${WORKING}/${sample}.unsorted.coverage.bed > \
     ${WORKING}/${sample}.coverage.bed
    done < ${list}

    ##Gate (20 second check; 5 minute report)
    #echo -e "\n"
    #GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_getcoverage" | wc -l)
    #GATEwait=0
    #until [[ $GATEcount == 0 ]]; do
    # sleep 20s
    # GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_getcoverage" | wc -l)
    # GATEwait=$[$GATEwait +1]
    # if [[ $GATEwait == 15 ]]; then
    #  echo "$(date): INCOMPLETE"
    #  echo "$(date): Waiting on ${GATEcount} coverage calculations to complete"
    #  GATEwait=0
    # fi
    #done
  fi

  #Merges output
  echo -e "\nCONSOLIDATING OUTPUT...\n"
  if [ ${mode} == "nucleotide" ]; then
    while read bam sample; do
      while read contig length; do
        cat ${WORKING}/${sample}.${contig}.coverage.bed >> ${WORKING}/${sample}.coverage.bed
      done < ${WORKING}/contig.list
    done < ${list}
  fi
  sed 1d ${list} > ${WORKING}/appendedsamples.tmp
  first=$( head -n1 ${list} | awk '{ print $2 }' )
  echo -e "Chr\nStart\nStop\n${first}" > ${WORKING}/${ID}.results.header.txt
  cp ${WORKING}/${first}.coverage.bed ${WORKING}/${ID}.${mode}.cov_matrix.bed
  while read bam sample; do
    echo "...appending ${sample}..."
    awk '{ print $4 }' ${WORKING}/${sample}.coverage.bed > ${WORKING}/${sample}.cov_vals.txt
    paste ${WORKING}/${ID}.${mode}.cov_matrix.bed ${WORKING}/${sample}.cov_vals.txt > \
    ${WORKING}/matrix.build.tmp
    mv ${WORKING}/matrix.build.tmp ${WORKING}/${ID}.${mode}.cov_matrix.bed
    echo "${sample}" >> ${WORKING}/${ID}.results.header.txt
  done < ${WORKING}/appendedsamples.tmp

  #Cleans and exits
  echo -e "\nCLEANING AND EXITING...\n"
  paste -s ${WORKING}/${ID}.results.header.txt > ${WORKING}/${ID}.results.header.reformatted.txt
  cat ${WORKING}/${ID}.results.header.reformatted.txt ${WORKING}/${ID}.${mode}.cov_matrix.bed > \
  ${OUTDIR}/${ID}.${mode}.cov_matrix.bed
  while read bam sample; do
    cp ${WORKING}/${sample}.coverage.bed ${OUTDIR}/raw_coverages/${sample}.coverage.bed
  done < ${list}
  sleep 10
  rm -rf ${WORKING}

  echo -e "COMPLETE\nFinal output written to:"
  echo "${OUTDIR}/${ID}.${mode}.cov_matrix.bed"

else
  echo -e "\n\nBinwise Coverage Calculator\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
  echo "Usage:"
  echo "  binCov.sh [samples.list] [ref.fa.dict] [nucleotide/physical] [binsize] [group_ID] [OUTDIR]"
  echo ""
fi


