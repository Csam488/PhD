#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#rpc_single.sh

set -e

usage(){
cat <<EOF
Usage: rpc_single.sh [-b] [-n N_MAD] (-d DIST | -m METRICS) SAMPLE BAM

SAMPLE, BAM, and one of DIST and METRICS are required arguments. 

If DIST is not provided, median insert size and MAD will be parsed from METRICS
to compute clustering distance. Distance is computed as median + N_MAD * MAD.

SAMPLE      Sample ID
BAM         Bam file (coordinate sorted)
-d DIST     Clustering distance (generally median insert + 7 MAD)
-m METRICS  Metrics file output by Picard's InsertSizeMetrics
-n N_MAD    Number of MAD added to median insert to compute minimum
            clustering distance [7]
-b          Convert rpc cluster files to bedpe format [FALSE]
            (Requires compress_clusters.py script from pycluster repo)

-h          Print this message
EOF
}

dist=""
metrics=""
bedpe=false
n_mad=7
while getopts ":d:m:bn:h" opt; do
  case "$opt" in
    d)
      dist=$OPTARG
      ;;
    m)
      metrics=$OPTARG
      ;;
    b)
      bedpe=true
      ;;
    n)
      n_mad=$OPTARG
      ;;
    h)
      usage
      exit 0
      ;;
  esac
done
shift $(( OPTIND - 1))

sample=$1
bam=$2

# If distance wasn't specified, calculate it from metrics file
if [[ -z $dist ]]; then
  if [[ -z $metrics ]]; then
    echo "ERROR: Must specify either distance or metrics file"
    usage
    exit 1
  fi

  read median mad <<<$(head -n8 $metrics | tail -n1 | cut -f1-2)
  dist=$(($median + $n_mad * $mad))
fi

# Extract discordant reads from bam
# Equivalent to `samtools view -F 3342`
sambamba view -f bam -F "not (proper_pair or unmapped or mate_is_unmapped or 
                              secondary_alignment or supplementary or 
                              duplicate)" $bam > ${sample}.disc.bam

# Convert discordant reads to bamstat format
# (Filters pairs where mapq=0 on both sides)
disc_to_rpc.py ${sample}.disc.bam ${sample}

for sv in del dup inv tloc; do
  python rpc.py \
    -d $dist \
    -x /scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/b37.lumpy.exclude.4-13.bed.gz \
    ${sample}.${sv}.pairs.txt \
    ${sample}.${sv}.clusters.txt
done

if [[ $bedpe ]]; then
  for sv in del dup inv tloc; do
    python compress_clusters.py \
      --no-split \
      ${sample}.${sv}.clusters.txt \
      ${sample}.${sv}.rpc.bedpe \
      ${sample}_${sv}
  done
fi
