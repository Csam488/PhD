#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Example parameters file
#Copy this format exactly for your own parameters file

###Run preferences###
#cohort ID (string)
export COHORT_ID=NZ 
#output directory - only final outputs will be written here
export OUTDIR=/scale_wlg_persistent/filesets/project/uoa02608/NZ_OUT
#working directory - all temp files will be written here. May temporarily need 
#lots of storage, so best to avoid writing to /tmp on restricted systems
export WRKDIR=${OUTDIR}/TMPFILES
#sex group to assign "other" sex for depth-based calling. Must be either "MALE" 
#or "FEMALE"
export other_assign=FEMALE
#underscore skip; determines how many underscores in clustering read ID are 
#assumed to be a part of sample ID; e.g. format like "STUDY_SAMPLE" 
#should have uscore_skip=1, etc
export uscore_skip=1
#Boolean (all caps) to indicate whether working directory should be kept once 
#pipeline is complete. Individual modules may be rerun with identical parameters
#file if working directory is saved.
export KEEP_TMP=FALSE
#boolean (all caps) to indicate if bamstat on all samples has already been run. 
#If set as TRUE, bamstat_paths MUST BE CORRECTLY SPECIFIED otherwise the 
#whole pipeline will fail
export pre_bamstat=FALSE 
#boolean (all caps) to indicate if you want to manually disable CNV genotyping. 
#Any value other than "TRUE" will revert to default operation (genotyping used 
#if cohort size >= ${min_geno}, set in mod. 6)
export GENOTYPE_OVERRIDE=FALSE
#maximum MAF/VAF before variant is considered artifactual and/or reference variant
export polyArt_filter=0.5

###Full paths for reference files & executables###
#path to Holmes git repo
export liWGS_SV=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes
#reference fasta
export REF=/scale_wlg_persistent/filesets/project/uoa02608/References/human_g2k_v37_ERCC92.fasta
#reference dictionary, restricted to chromosomes where calls should be made
export DICT=/scale_wlg_persistent/filesets/project/uoa02608/References/human_g2k_v37_ERCC92.dict
#classifier git repo
export CLASSIFIER_DIR=${liWGS_SV}/classifier
#pycluster git repo
export PYCLUSTER_DIR=${liWGS_SV}/pycluster
#picard.jar executable
export PICARD=/scale_wlg_persistent/filesets/opt_nesi/mahuika/picard/2.1.0/picard.jar
#sambamba executable
export sambamba=/scale_wlg_persistent/filesets/opt_nesi/mahuika/Sambamba/0.6.5-gimkl-2017a/bin/sambamba
#blacklist file, restrict CNV calling on >30% coverage of features in list 
export CNV_BLACKLIST=/scale_wlg_persistent/filesets/project/uoa02608/consensusBlacklist.bed
#file containing paths to pre-run bamstat directories. First column: ID, 
#second column: full path to bamstat directory; tab-delimited. Will be 
#ignored unless pre_bamstat="TRUE"
export bamstat_paths=${OUTDIR}/bamstat.list
#path to antibody parts annotation file (available from UCSC), 
#used for exclusion of putative complex site
export abParts=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/abParts.bed
#bed file corresponding to N-masked regions of reference genome
export NMASK=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/hg19_N_masked.bed
#UCSC refFlat (genes)
export refFlat=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/refFlat.bed
#Gencode GTF
export GTF=/scale_wlg_persistent/filesets/project/uoa02608/modules/Holmes/data/gencode.v29.annotation.gtf
#TEMP dir
export TMPDIR=/scale_wlg_nobackup/filesets/nobackup/uoa02608

#Update user paths
export PATH=${PATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}
export PYTHONPATH=${PYTHONPATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}