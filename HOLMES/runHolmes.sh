#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Submit this script to run entire pipeline

#Read input
samples_list=$1
params=$2
part=$3 #A=mods 1-5, B=mods 6-end, F=full pipeline
#Ensure correct version of gcc for anaconda dependencies

module load gcc/6.3.0
module load Sambamba/0.6.5-gimkl-2017a
module load BEDTools/2.26.0-gimkl-2017a
module load BamTools/2.4.1-gimkl-2017a
module load SAMtools/1.8-gimkl-2017a
module load Python/2.7.14-gimkl-2017a
module load Perl/5.24.1-gimkl-2017a
module load GLib/2.52.0-gimkl-2017a
module load BLAT/3.5-gimkl-2017a
#module load BioConductor/3.5-gimkl-2017a-R-3.4.2
#module load R/3.5.0-gimkl-2017a
module load R/3.6.1-gimkl-2018b
#################################NEED ANACONDA?#################################

echo -e "STATUS [$(date)]: Loading parameters..."
. ${params}
#Source params file
#####FIX INSTALL DIRECTORY#####
export R_LIBS=${liWGS_SV}/R
export R_LIBS_USER=${liWGS_SV}/R

#Skip modules 1-5 if part=B
if [ ${part} == "A" ] || [ ${part} == "F" ]; then
  #Set up output directory tree
  echo -e "STATUS [$(date)]: Creating directory trees..."
  if ! [ -e ${OUTDIR} ]; then
    mkdir ${OUTDIR}
  fi
  if ! [ -e ${TMPDIR} ]; then
    mkdir ${TMPDIR}
  fi
  cp ${samples_list} ${OUTDIR}/samples.list
  samples_list=$(echo -e "${OUTDIR}/samples.list") ###Copy sample list to prevent direcory errors later
  if ! [ -e ${liWGS_SV}/R ]; then
    mkdir ${liWGS_SV}/R
  fi
  mkdir ${OUTDIR}/checkpoints
  mkdir ${OUTDIR}/QC
  mkdir ${OUTDIR}/QC/cohort
  mkdir ${OUTDIR}/QC/sample
  while read ID bam sex; do
    mkdir ${OUTDIR}/QC/sample/${ID}
  done < ${samples_list}
  mkdir ${OUTDIR}/data
  mkdir ${OUTDIR}/data/seqDepth
  mkdir ${OUTDIR}/data/clusters
  mkdir ${OUTDIR}/SV_calls
  mkdir ${OUTDIR}/logs
  mkdir ${OUTDIR}/plots
  cp ${params} ${OUTDIR}/${COHORT_ID}.run_parameters.sh

  #Set up working directory tree
  if ! [ -e ${WRKDIR} ]; then
    mkdir ${WRKDIR}
  fi
  while read ID bam sex; do
    mkdir ${WRKDIR}/${ID}
  done < ${samples_list}
  #cp ${liWGS_SV}/scripts/SE_largeCNV/* ${WRKDIR}/

  #Write start date to temporary file
  echo $(date) > ${WRKDIR}/start.tmp

  #Make 2bit
  if ! [ -e ${REF}.2bit ]; then
    echo -e "STATUS [$(date)]: Making $REF.2bit..."
    /scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/faToTwoBit ${REF} ${REF}.2bit
  fi

  #Link & index all bams
  echo -e "STATUS [$(date)]: Indexing BAMs..."
  while read ID bam sex; do
    echo -e "STATUS [$(date)]: Indexing $ID..."
    if ! [ -e ${WRKDIR}/${ID}/${ID}.bam ]; then
      ln -s ${bam} ${WRKDIR}/${ID}/${ID}.bam
    fi
    if ! [ -e ${WRKDIR}/${ID}/${ID}.bam.bai ]; then
      sambamba index ${WRKDIR}/${ID}/${ID}.bam
    fi
    echo -e "STATUS [$(date)]: Indexed $ID..."
  done < ${samples_list}

  #Gate until indexing is complete; 20 sec check; 5 min report
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_index" | wc -l )
#  GATEwait=0
#  until [[ $GATEcount == 0 ]]; do
#    sleep 20s
#    GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_index" | wc -l )
#    GATEwait=$[${GATEwait} +1]
#    if [[ $GATEwait == 15 ]]; then
#      echo -e "STATUS [$(date)]: Gated at bam indexing..."
#      GATEwait=0
#    fi
#  done

  ##STAGE 1: modules 1, 2, and 3
  echo -e "STATUS [$(date)]: Beginning PHASE 1..."
  #Submit module 1 (QC)
  #bsub -q normal -sla miket_sc -o (out) ${OUTDIR}/logs/module1.log -e (error) ${OUTDIR}/logs/module1.log
  if ! [ -e ${OUTDIR}/checkpoints/module1.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 1..."
    bash "${liWGS_SV}/scripts/module1.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 1 found..."
  fi
  #Submit module 2 (physical depth analysis)
  if ! [ -e ${OUTDIR}/checkpoints/module2.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 2..."
    bash "${liWGS_SV}/scripts/module2.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 2 found..."
  fi
  #Submit module 3 (per-sample clustering)
  if ! [ -e ${OUTDIR}/checkpoints/module3.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 3..."
    bash "${liWGS_SV}/scripts/module3.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 3 found..."
  fi

  #Gate until modules 1 & 2 complete; 20 sec check; 5 min report
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_MODULE_2" | wc -l )
#  GATEwait=0
#  until [[ $GATEcount == 0 ]]; do
#    sleep 20s
#    GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_MODULE_2" | wc -l )
#    GATEwait=$[${GATEwait} +1]
#    if [[ $GATEwait == 15 ]]; then
#      echo -e "STATUS [$(date)]: Gated at PHASE 1a..."
#      GATEwait=0
#    fi
#  done

  ##STAGE 2a: module 4
  echo -e "STATUS [$(date)]: PHASE 1a complete; Beginning PHASE 2a..."
  #Submit module 4 (physical depth CNV calling)
  if ! [ -e ${OUTDIR}/checkpoints/module4.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 4..."
    bash "${liWGS_SV}/scripts/module4.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 4 found..."
  fi

  #Gate until module 3 complete; 20 sec check; 5 min report
  #GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_3" | wc -l )
  #GATEwait=0
  #until [[ $GATEcount == 0 ]]; do
  #  sleep 20s
  #  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_3" | wc -l )
  #  GATEwait=$[${GATEwait} +1]
  #  if [[ $GATEwait == 15 ]]; then
  #    echo -e "STATUS [$(date)]: Gated at PHASE 1b..."
  #    GATEwait=0
  #  fi
  #done

  ##STAGE 2b: module 5
  echo -e "STATUS [$(date)]: PHASE 1b complete; Beginning PHASE 2b..."
  #Submit module 5 (joint clustering)
  if ! [ -e ${OUTDIR}/checkpoints/module5.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 5..."
    bash "${liWGS_SV}/scripts/module5.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 5 found..."
  fi

  #Gate until modules 4 & 5 complete; 20 sec check; 5 min report
  #GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_4\|${COHORT_ID}_MODULE_5" | wc -l )
  #GATEwait=0
  #until [[ $GATEcount == 0 ]]; do
  #  sleep 20s
  #  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_4\|${COHORT_ID}_MODULE_5" | wc -l )
  #  GATEwait=$[${GATEwait} +1]
  #  if [[ $GATEwait == 15 ]]; then
  #    echo -e "STATUS [$(date)]: Gated at PHASE 2..."
  #    GATEwait=0
  #  fi
  #done
fi

#Skip modules 6-end if part A is optioned
if [ ${part} == "B" ] || [ ${part} == "F" ]; then
  ##STAGE 3: modules 6 and 7
  echo -e "STATUS [$(date)]: Beginning PHASE 3..."
  #Submit module 6 (consensus CNV merging)
  if ! [ -e ${OUTDIR}/checkpoints/module6.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 6..."
    bash "${liWGS_SV}/scripts/module6.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 6 found..."
  fi
  #Submit module 7 (complex SV categorization)
  if ! [ -e ${OUTDIR}/checkpoints/module7.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 7..."
    bash "${liWGS_SV}/scripts/module7.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 7 found..."
  fi

  #Gate until modules 6 & 7 completes; 20 sec check; 5 min report
  #GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_6\|${COHORT_ID}_MODULE_7" | wc -l )
  #GATEwait=0
  #until [[ $GATEcount == 0 ]]; do
  #  sleep 20s
  #  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_6\|${COHORT_ID}_MODULE_7" | wc -l )
  #  GATEwait=$[${GATEwait} +1]
  #  if [[ $GATEwait == 15 ]]; then
  #    echo -e "STATUS [$(date)]: Gated at PHASE 3..."
  #    GATEwait=0
  #  fi
  #done

  ##STAGE 4: module 8
  echo -e "STATUS [$(date)]: Beginning PHASE 4..."
  if ! [ -e ${OUTDIR}/checkpoints/module8.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 8..."
    bash "${liWGS_SV}/scripts/module8.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 8 found..."
  fi

  ##STAGE 5: module 9
  echo -e "STATUS [$(date)]: Beginning PHASE 5..."
  if ! [ -e ${OUTDIR}/checkpoints/module9.done ]; then
    echo -e "STATUS [$(date)]: ====================== Executing module 9..."
    bash "${liWGS_SV}/scripts/module9.sh" $samples_list $params
  else
    echo -e "STATUS [$(date)]: ====================== Checkpoint for module 9 found..."
  fi

  #Gate until module 9 completes; 20 sec check; 5 min report
  #GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_9" | wc -l )
  #GATEwait=0
  #until [[ $GATEcount == 0 ]]; do
  #  sleep 20s
  #  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_9" | wc -l )
  #  GATEwait=$[${GATEwait} +1]
  #  if [[ $GATEwait == 15 ]]; then
  #    echo -e "STATUS [$(date)]: Gated at PHASE 5..."
  #    GATEwait=0
  #  fi
  #done

  #Move all relevant files to OUTDIR
  cp ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${OUTDIR}/data/seqDepth/
  cp ${WRKDIR}/final_variants/* ${OUTDIR}/SV_calls
  cp ${WRKDIR}/raw_clusters/* ${OUTDIR}/data/clusters
  mkdir ${OUTDIR}/SV_calls/annotations
  cp ${WRKDIR}/annotations/*_gene_anno.bed ${OUTDIR}/SV_calls/annotations/

  #Summarize run outcomes
  echo -e "STATUS [$(date)]: Calculating run summary metrics..."
  bash ${liWGS_SV}/scripts/Holmes_summary.sh ${samples_list} ${params}

  #Remove working directory unless specified otherwise
  if ! [ ${KEEP_TMP}=="TRUE" ]; then
    rm -rf ${WRKDIR}
  fi

  #Print final completion notice
  echo -e "-----\nSTATUS [$(date)]: Holmes liWGS-SV discovery on ${COHORT_ID} complete\n-----\n\n"
fi
