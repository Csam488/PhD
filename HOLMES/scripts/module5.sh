#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Module 5 (Joint Clustering & Classification)

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Prepare classifier files
mkdir ${WRKDIR}/classifier
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/bamstat"
done < ${samples_list} > ${WRKDIR}/classifier/classifier.samples.list
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.coverage.bed.gz"
done < ${samples_list} > ${WRKDIR}/classifier/classifier.icov.list

#Run classifier
echo -e "STATUS [$(date)]: Running Classifier..."
#cd ${WRKDIR}/classifier
bash ${CLASSIFIER_DIR}/run_classify.sh -m ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ${WRKDIR}/classifier/classifier.samples.list ${COHORT_ID} ${WRKDIR}/classifier/classifier.icov.list

#Gate until complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo "$(date): INCOMPLETE"
#    echo "$(date): Waiting on ${GATEcount} jobs to complete"
#    GATEwait=0
#  fi
#done

#Apply cluster fix
echo -e "STATUS [$(date)]: Applying Cluster Fix..."
mkdir ${WRKDIR}/classifier/clusterfix
mkdir ${WRKDIR}/classifier/clusterfix/delsplit
for contig in $( seq 1 22 ) X Y; do
  echo "STATUS [$(date)]: For Ch ${contig}"
  awk -v OFS="\t" -v contig=${contig} '{ if ($4==contig || $4=="") print $0 }' ${WRKDIR}/classifier/${COHORT_ID}_deletion.clusters | cat -s | sed '1{/^$/d}' > ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.clusters
  bash ${CLASSIFIER_DIR}/cleanClusters_patch.sh ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.clusters ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.patched.clusters ${uscore_skip} TRUE
done
for class in insertion inversion transloc; do
  echo "STATUS [$(date)]: For ${class}"
  bash ${CLASSIFIER_DIR}/cleanClusters_patch.sh ${WRKDIR}/classifier/${COHORT_ID}_${class}.clusters ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${uscore_skip} FALSE
done

#Gate until complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_clusterPatch" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_clusterPatch" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo "$(date): INCOMPLETE"
#    echo "$(date): Waiting on ${GATEcount} jobs to complete"
#    GATEwait=0
#  fi
#done

#Concatenate split deletion patched clusters
echo -e "STATUS [$(date)]: Concatenating Split Clusters..."
for contig in $( seq 1 22 ) X Y; do
  echo "STATUS [$(date)]: For Ch ${contig}"
  cat ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.patched.clusters >> ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_deletion.patched.clusters
done

#Complete classifier on fixed clusters
echo -e "STATUS [$(date)]: Completing Classifier..."
for class in deletion insertion inversion transloc; do
  echo "STATUS [$(date)]: For ${class}"
  python ${CLASSIFIER_DIR}/rpc_classify.py ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.events.bedpe ${WRKDIR}/classifier/${COHORT_ID}_boot.list ${COHORT_ID}_${class} --cluster-bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.bkpts.bedpe
done

#Gate until complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier_r2" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier_r2" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo "$(date): INCOMPLETE"
#    echo "$(date): Waiting on ${GATEcount} jobs to complete"
#    GATEwait=0
#  fi
#done

#Reclassify fixed clusters
echo -e "STATUS [$(date)]: Reclassifying Fixed Clusters..."
mkdir ${WRKDIR}/classifier/clusterfix/newCoords
for class in deletion insertion inversion transloc; do
  bash ${CLASSIFIER_DIR}/reclassify_output.sh ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.events.bedpe ${class} ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${WRKDIR}/classifier/clusterfix/newCoords FALSE
done

#Gate until complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_reclassify" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_reclassify" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo "$(date): INCOMPLETE"
#    echo "$(date): Waiting on ${GATEcount} jobs to complete"
#    GATEwait=0
#  fi
#done
echo -e "$(date)" > ${OUTDIR}/checkpoints/module5.done