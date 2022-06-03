#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Module 7 (Complex SV Categorization)

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Dynamically determine min complex size
minCXsize=$( fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | cut -f9 | sort -nrk1,1 | head -n1 )

#Make cnMOPS list
if [ -e ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list ]; then
  rm ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list
fi
if [ -e ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list ]; then
  rm ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list
fi
while read ID bam sex; do
	echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed" >> ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list
	echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed" >> ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list
done < ${samples_list}

#Submit complex linking script
echo -e "STATUS [$(date)]: Performing complex linking"
bash ${liWGS_SV}/scripts/complexLinking.sh ${samples_list} ${params}

#Submit inversion categorization
echo -e "STATUS [$(date)]: Performing inversion categorization"
bash ${liWGS_SV}/scripts/classify_inversion.sh ${WRKDIR}/classifier/clusterfix/newCoords/inversion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_inversion.patched.clusters ${minCXsize} ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list ${WRKDIR}/classifier/clusterfix/newCoords/

#Submit translocation categorization
echo -e "STATUS [$(date)]: Performing translocation categorization"
bash ${liWGS_SV}/scripts/classify_translocation.sh ${WRKDIR}/classifier/clusterfix/newCoords/transloc.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_transloc.patched.clusters ${WRKDIR}/classifier/clusterfix/newCoords/

#Gate until complex linking complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cxLINK" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cxLINK" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo -e "STATUS [$(date)]: Waiting on complex linking script..."
#    GATEwait=0
#  fi
#done

#Determine clustering distance (max of all samples)
clustdist=$( cut -f3 ${WRKDIR}/classifier/${COHORT_ID}_boot.list | sort -nrk1,1 | head -n1 )

#Parse complex linked output
echo -e "STATUS [$(date)]: Parsing output"
if [ -e ${WRKDIR}/events.list ]; then
  rm ${WRKDIR}/events.list
fi
if [ -e ${WRKDIR}/clusters.list ]; then
  rm ${WRKDIR}/clusters.list
fi
for type in deletion insertion inversion transloc; do
  echo -e "${type}\t${WRKDIR}/classifier/clusterfix/newCoords/${type}.events.reclassified.bedpe" >> ${WRKDIR}/events.list
  echo -e "${type}\t${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${type}.patched.clusters" >> ${WRKDIR}/clusters.list
done
bash ${liWGS_SV}/scripts/parseComplex.sh ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list ${WRKDIR}/events.list ${WRKDIR}/clusters.list ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list ${WRKDIR}/classifier/clusterfix/newCoords/ ${clustdist} $( cat ${samples_list} | wc -l ) ${params}

#Gate until complex linking complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classification" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classification" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo -e "STATUS [$(date)]: Waiting on BCA categorization..."
#    GATEwait=0
#  fi
#done
echo -e "$(date)" > ${OUTDIR}/checkpoints/module7.done