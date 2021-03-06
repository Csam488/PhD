#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Module 9 (Variant annotation)

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Make output directory
mkdir ${WRKDIR}/annotations

#Submit deletion annotation
echo -e "STATUS [$(date)]: Submitting Deletion Annotation"
awk -v OFS="\t" '{ print $1, $2, $3, "DEL", $4 }' ${WRKDIR}/final_variants/${COHORT_ID}.deletion.bed | fgrep -v "#" > ${WRKDIR}/annotations/deletion_preAnno.bed
bash ${liWGS_SV}/scripts/annotate_SVintervals.sh ${WRKDIR}/annotations/deletion_preAnno.bed DEL ${WRKDIR}/annotations/deletion_gene_anno.bed ${params}

#Submit duplication annotation
echo -e "STATUS [$(date)]: Submitting Duplication Annotation"
awk -v OFS="\t" '{ print $1, $2, $3, "DUP", $4 }' ${WRKDIR}/final_variants/${COHORT_ID}.duplication.bed | fgrep -v "#" > ${WRKDIR}/annotations/duplication_preAnno.bed
bash ${liWGS_SV}/scripts/annotate_SVintervals.sh ${WRKDIR}/annotations/duplication_preAnno.bed DUP ${WRKDIR}/annotations/duplication_gene_anno.bed ${params}

#Submit inversion annotation (averages cluster coordinates for breakpoints in BED)
echo -e "STATUS [$(date)]: Submitting Inversion Annotation"
awk -v OFS="\t" '{ printf "%s\t%.0f\t%.0f\t%s\t%s\n", $1, ($2+$5)/2, ($3+$6)/2, "INV", $7 }' ${WRKDIR}/final_variants/${COHORT_ID}.inversion.bedpe | fgrep -v "#" > ${WRKDIR}/annotations/inversion_preAnno.bed
bash ${liWGS_SV}/scripts/annotate_SVintervals.sh ${WRKDIR}/annotations/inversion_preAnno.bed INV ${WRKDIR}/annotations/inversion_gene_anno.bed ${params}

#Submit insertion source annotation
echo -e "STATUS [$(date)]: Submitting Insertion Source Annotation"
awk -v OFS="\t" '{ print $1, $2, $3, "INS_SOURCE", $7 }' ${WRKDIR}/final_variants/${COHORT_ID}.insertion.bedpe | fgrep -v "#" > ${WRKDIR}/annotations/insertionSource_preAnno.bed
bash ${liWGS_SV}/scripts/annotate_SVintervals.sh ${WRKDIR}/annotations/insertionSource_preAnno.bed INS_SOURCE ${WRKDIR}/annotations/insertionSource_gene_anno.bed ${params}

#Submit insertion sink annotation
echo -e "STATUS [$(date)]: Submitting Insertion Sink Annotation"
awk -v OFS="\t" '{ print $4, $5, $6, "INS_SINK", $7 }' ${WRKDIR}/final_variants/${COHORT_ID}.insertion.bedpe | fgrep -v "#" | sed '1d' > ${WRKDIR}/annotations/insertionSink_preAnno.bed
bash ${liWGS_SV}/scripts/annotate_SVintervals.sh ${WRKDIR}/annotations/insertionSink_preAnno.bed INS_SINK ${WRKDIR}/annotations/insertionSink_gene_anno.bed ${params}

#Submit retrotransposon check
echo -e "STATUS [$(date)]: Checking for Retrotransposons"
for class in deletion insertion inversion transloc; do
  echo -e "${class}\t${WRKDIR}/classifier/clusterfix/newCoords/${class}.events.reclassified.bedpe" >> ${WRKDIR}/events.list
  echo -e "${class}\t${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters" >> ${WRKDIR}/clusters.list
done
bash ${liWGS_SV}/scripts/get_retrotransposons.sh ${WRKDIR}/events.list ${WRKDIR}/clusters.list ${WRKDIR}/classifier/clusterfix/newCoords/

#Gate until all annotations complete; 20 sec check; 5 min report
#GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_annotation" | wc -l )
#GATEwait=0
#until [[ $GATEcount == 0 ]]; do
#  sleep 20s
#  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_annotation" | wc -l )
#  GATEwait=$[${GATEwait} +1]
#  if [[ $GATEwait == 15 ]]; then
#    echo -e "STATUS [$(date)]: Waiting on annotations..."
#    GATEwait=0
#  fi
#done
echo -e "$(date)" > ${OUTDIR}/checkpoints/module9.done