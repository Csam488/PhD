#!/bin/bash

#################################
# Converter and annotation      #
# script for use with GRIDSS    #
#################################

input=$1
output=$2
freq=$3
max_occ=$4
annotation_dir=$5

min_samples=1000 #min samples for DGV freq

####################################################
echo -e "STATUS [$(date)]: Starting..."
####################################################

module load Python/2.7.14-gimkl-2017a
module load BEDTools/2.26.0-gimkl-2017a
module load Tcl/8.6.6-gimkl-2017a

function join_by { local IFS="$1"; shift; echo "$*"; }

if ! [ -e ${output} ]; then
	mkdir ${output}
fi
if ! [ -e ${output}/Annotat ]; then
	mkdir ${output}/Annotat
fi

####################################################
echo -e "STATUS [$(date)]: Comparing to Internal Dataset..."
####################################################
grep "DUP" ${input} | awk -v T='\t' -F '\t' '{if ($3 < $4 ) print $2T$3T$4T$8T$9T$10T$11; else print $2T$4T$3T$8T$9T$10T$11}' > ${output}/DUP.bed
grep "DEL" ${input} | awk -v T='\t' -F '\t' '{if ($3 < $4 ) print $2T$3T$4T$8T$9T$10T$11; else print $2T$4T$3T$8T$9T$10T$11}' > ${output}/DEL.bed
grep "INV" ${input} | cut -f 1,2,3,4,6,7,8,9,10,11 | awk -v T='\t' -F '\t' '{N=$3; M=$3; for (i=4; i<=6; i++) if ($i <N) N= $i; for (i=4; i<=6; i++) if ($i >M) M= $i; print $2T N T M T$7T$8T$9T$10}' > ${output}/INV.bed
grep "TRLX" ${input} | awk -v T='\t' -F '\t' '{if ($3 < $4 ) print $2T$3T$4T$5T$6T$7T$8T$9T$10T$11; else print $2T$4T$3T$5T$6T$7T$8T$9T$10T$11}' > ${output}/TRLX.bed
awk -v T='\t' -F '\t' '{if ($5 < $6 ) print $1T$2T$3T$4T$5T$6T$7T$8T$9T$10; else print $1T$2T$3T$4T$6T$5T$7T$8T$9T$10}' ${output}/TRLX.bed > ${output}/TRLX.tmp.bed
awk -v T='\t' -v X='X' -F '\t' '{if ($1 < $4 && $1 != X) print $1T$2T$3T$4T$5T$6T$7T$8T$9T$10; else print $4T$5T$6T$1T$2T$3T$7T$8T$9T$10}' ${output}/TRLX.tmp.bed | grep -v "ORPHAN" > ${output}/TRLX.bed
#awk -v T='\t' -v X='X' -F '\t' '{if ($1 < $4 && $1 != X) print $1T$2T$3T$4T$5T$6T$7T$8T$9T$10; else print $4T$5T$6T$1T$2T$3T$7T$8T$9T$10}' ${output}/TRLX.tmp.bed | grep "ORPHAN" > ${output}/ORPHAN_TRLX.bed
rm ${output}/TRLX.tmp.bed

if [ -e ${output}/ALL.unsort.unann.bed ]; then
	rm ${output}/ALL.unsort.unann.bed
fi
echo "0" > ${output}/0.tmp
while read Chr Start Stop Other; do
	echo -e "${Chr}\t${Start}\t${Stop}" > ${output}/Temp.bed
	echo -e "${Chr}\t${Start}\t${Stop}\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${annotation_dir}/Data/DUP_38.bed -b ${output}/Temp.bed | sort -k4,4n | cut -f 4 | cat - ${output}/0.tmp | head -n 1)" >> ${output}/ALL.unsort.unann.bed
done < ${output}/DUP.bed
while read Chr Start Stop Other; do
	echo -e "${Chr}\t${Start}\t${Stop}" > ${output}/Temp.bed
	echo -e "${Chr}\t${Start}\t${Stop}\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${annotation_dir}/Data/DEL_38.bed -b ${output}/Temp.bed | sort -k4,4n | cut -f 4 | cat - ${output}/0.tmp | head -n 1)" >> ${output}/ALL.unsort.unann.bed
done < ${output}/DEL.bed
while read Chr Start Stop Other; do
	echo -e "${Chr}\t${Start}\t${Stop}" > ${output}/Temp.bed
	echo -e "${Chr}\t${Start}\t${Stop}\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${annotation_dir}/Data/INV_38.bed -b ${output}/Temp.bed | sort -k4,4n | cut -f 4 | cat - ${output}/0.tmp | head -n 1)" >> ${output}/ALL.unsort.unann.bed
done < ${output}/INV.bed
while read Chr Start Stop Chr2 Start2 Stop2 Name Other; do
	echo -e "${Chr}\t$(( ${Start} - 500))\t$(( ${Stop} + 500))" > ${output}/Temp.bed
	bedtools intersect -f 0.7 -F 0.7 -u -a ${annotation_dir}/Data/TRLX_38.bedpe -b ${output}/Temp.bed | awk -v T='\t' -F '\t' '{print $4T$5T$6T$1T$2T$3T$7}' > ${output}/Temp.Match.bed
	echo -e "${Chr2}\t$(( ${Start2} - 500))\t$(( ${Stop2} + 500))" > ${output}/Temp.bed
	echo -e "${Chr}\t${Start}\t${Stop}\t${Name}_A\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${output}/Temp.Match.bed -b ${output}/Temp.bed | sort -k4,4n | cat - ${output}/0.tmp | cut -f 7 | head -n 1)" >> ${output}/ALL.unsort.unann.bed
	echo -e "${Chr2}\t${Start2}\t${Stop2}\t${Name}_B\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${output}/Temp.Match.bed -b ${output}/Temp.bed | sort -k4,4n | cat - ${output}/0.tmp | cut -f 7 | head -n 1)" >> ${output}/ALL.unsort.unann.bed
done < ${output}/TRLX.bed
#while read Chr Start Stop Chr2 Start2 Stop2 Name Other; do
#	echo -e "${Chr}\t$(( ${Start} - 500))\t$(( ${Stop} + 500))" > ${output}/Temp.bed
#	bedtools intersect -f 0.7 -F 0.7 -u -a ${annotation_dir}/Data/ORPHAN_TRLX.bedpe -b ${output}/Temp.bed | awk -v T='\t' -F '\t' '{print $4T$5T$6T$1T$2T$3T$7}' > ${output}/Temp.Match.bed
#	echo -e "${Chr2}\t$(( ${Start2} - 500))\t$(( ${Stop2} + 500))" > ${output}/Temp.bed
#	echo -e "${Chr}\t${Start}\t${Stop}\t${Name}_A\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${output}/Temp.Match.bed -b ${output}/Temp.bed | sort -k4,4n | cat - ${output}/0.tmp | cut -f 4 | head -n 1)" >> ${output}/ALL.unsort.unann.bed
#	echo -e "${Chr2}\t${Start2}\t${Stop2}\t${Name}_B\t${Other}\t$(bedtools intersect -f 0.7 -F 0.7 -u -a ${output}/Temp.Match.bed -b ${output}/Temp.bed | sort -k4,4n | cat - ${output}/0.tmp | cut -f 4 | head -n 1)" >> ${output}/ALL.unsort.unann.bed
#done < ${output}/ORPHAN_TRLX.bed
rm ${output}/Temp.bed ${output}/0.tmp
awk -v OFS="\t" -v T=$max_occ '{if ( $8 <= T ) print}' ${output}/ALL.unsort.unann.bed > ${output}/ALL.unsort.filt.unann.bed

echo -e "STATUS [$(date)]:\t\t$(( $(wc -l ${output}/ALL.unsort.unann.bed | cut -d ' ' -f 1) - $(wc -l ${output}/ALL.unsort.filt.unann.bed | cut -d ' ' -f 1) )) records removed as they exceed an incidence of ${max_occ} including $(( $(grep "TRLX" ${output}/ALL.unsort.unann.bed | wc -l) - $(grep "TRLX" ${output}/ALL.unsort.filt.unann.bed | wc -l) )) TRLX records"










####################################################
echo -e "STATUS [$(date)]: Making master BED file..."
####################################################
if [ -e ${output}/ALL.sort.unann.bed ]; then
	rm ${output}/ALL.sort.unann.bed
fi
for chr in $( seq 1 22 ) X Y; do
	awk -v OFS="\t" -v CHR=$chr '{if ($1==CHR) print}' ${output}/ALL.unsort.filt.unann.bed |sort -k1,1n -k2,2n -k3,3n >> ${output}/ALL.sort.unann.bed
done
rm ${output}/ALL.unsort.unann.bed
rm ${output}/ALL.unsort.filt.unann.bed


####################################################
echo -e "STATUS [$(date)]: Running AnnotSV..."
####################################################

export ANNOTSV=/scale_wlg_persistent/filesets/project/uoa02608/modules/AnnotSV_2.1
$ANNOTSV/bin/AnnotSV -genomeBuild GRCh38 -SVinputFile ${output}/ALL.sort.unann.bed -outputdir ${output}/Annotat > ${output}/AnnotatSV.log
sub="$(fgrep -v "#" ${output}/ALL.sort.unann.bed|wc -l|head -n 1)"
ann="$(awk -v OFS="\t" -v type="full" '{if ($11==type) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv|wc -l |head -n 1)"
####################################################
echo -e "STATUS [$(date)]: Compiling annotation..."
####################################################

if [ ${sub} != ${ann} ]; then
	echo -e "ERROR:: annotated records does not equal submitted records: ${sub} , ${ann}"
	exit
else
	#echo -e "#Chr\tStart\tStop\tType\tName\tScore\tFilter\tInternal_Count\tDGV_gain_IDs\tDGV_gain_freq\tDGV_loss_IDs\tDGV_loss_freq\tsp_DGV_gain_IDs\tsp_DGV_gain_freq\tsp_DGV_loss_IDs\tsp_DGV_loss_freq\tDDD_SVs\tDDD_gain_freq\tDDD_loss_freq\tsp_DDD_SVs\tsp_DDD_gain_freq\tsp_DDD_loss_freq\tkg_events\tkg_freq\tkg_maxfreq\tsp_kg_events\tsp_kg_freq\tsp_kg_maxfreq\tGD_ID\tGD_N_HOM\tGD_AT\tGD_POPMAX\tGD_ID_others\tIMH_ID\tIMH_AF\tIMH_ID_others\tdbvar_events\tdbvar_IDs\tdbvar_status\tsp_dbvar_events\tsp_dbvar_IDs\tsp_dbvar_status\tTADcoordinates\tENCODE_expiriments\tGC_left\tGC_right\trepeat_loc_left\trepeat_type_left\trepeat_loc_right\trepeat_type_right\tgene\tnm\tCDS_len\ttx_len\tloc\tint_st\tint_end\tpromoter\tHI_GC\tTriS_GC\tDDD_stat\tDDD_mode\tDDD_disease\tDDD_pmid\tDDD_HI\tExAC_synZ\tExAC_misZ\tExAC_pLI\tExAC_delz\tExAC_dupz\tExAC_cnvz\tmorbidGene\tmorbidGene_candidate\tMIM\tPhenotypes\tInheritance" > ${output}/ALL.sort.ann.bed
	echo -e "#Chr\tStart\tStop\tType\tName\tScore\tFilter\tInternal_Count\tDGV_gain_IDs\tDGV_gain_freq\tDGV_loss_IDs\tDGV_loss_freq\tsp_DGV_gain_IDs\tsp_DGV_gain_freq\tsp_DGV_loss_IDs\tsp_DGV_loss_freq\tkg_events\tkg_freq\tkg_maxfreq\tsp_kg_events\tsp_kg_freq\tsp_kg_maxfreq\tIMH_ID\tIMH_AF\tIMH_ID_others\tdbvar_events\tdbvar_IDs\tdbvar_status\tsp_dbvar_events\tsp_dbvar_IDs\tsp_dbvar_status\tTADcoordinates\tENCODE_expiriments\tGC_left\tGC_right\trepeat_loc_left\trepeat_type_left\trepeat_loc_right\trepeat_type_right\tgene\tCDS_len\ttx_len\tloc\tint_st\tint_end\tpromoter\tHI_GC\tTriS_GC\tDDD_stat\tDDD_mode\tDDD_disease\tDDD_pmid\tDDD_HI\tExAC_synZ\tExAC_misZ\tExAC_pLI\tExAC_delz\tExAC_dupz\tExAC_cnvz\tmorbidGene\tmorbidGene_candidate\tMIM\tPhenotypes\tInheritance" > ${output}/ALL.sort.ann.bed
	while read Chr Start Stop Name Score Filter Type Freq ; do
		awk -v OFS="\t" -v type="split" -v name=${Name} '{if ($11==type && $6==name) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/\t\t/\t-\t/g'| sed 's/\t\t/\t-\t/g'| sed 's/ /_/g'| sed 's/\t-1/\t-/g' | sort -u > ${output}/${Name}.split.tmp
		#awk -v OFS="\t" -v type="full" -v name=${Name} '{if ($11==type && $6==name) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/\t\t/\t-\t/g'| sed 's/\t\t/\t-\t/g'| sed 's/ /_/g'| sed 's/\t-1/\t-/g' | sort -u > ${output}/${Name}.full.tmp ###HEAD HERE IF ISSUE WITH ANNOTAT CONTINUES
		awk -v OFS="\t" -v type="full" -v name=${Name} '{if ($11==type && $6==name) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv | head -n 1| sed 's/\t\t/\t-\t/g'| sed 's/\t\t/\t-\t/g'| sed 's/ /_/g'| sed 's/\t-1/\t-/g' | sort -u > ${output}/${Name}.full.tmp ###HEAD HERE IF ISSUE WITH ANNOTAT CONTINUES
		if [ $(wc -l < ${output}/${Name}.split.tmp | head -n 1) == "0" ]; then
			echo -e "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-" > ${output}/${Name}.split.tmp
		fi
		DGV_gain_IDs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_IDs | cut -d: -f1) ${output}/${Name}.full.tmp)
		DGV_gain_freq=$( cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_n_samples_tested | cut -d: -f1),$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | awk -v min=${min_samples} '{if ( $1>=min ) print $2 ; else print "Fewer_than_"min"_samples"}' | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		DGV_loss_IDs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_IDs | cut -d: -f1) ${output}/${Name}.full.tmp)
		DGV_loss_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_n_samples_tested | cut -d: -f1),$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | awk -v OFS="\t"  -v min=${min_samples} '{if ( $1>=min ) print $2 ; else print "Fewer_than_"min"_samples"}' | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		#DDD_SVs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_SV | cut -d: -f1) ${output}/${Name}.full.tmp)
		#DDD_gain_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DUP_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		#DDD_loss_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DEL_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		kg_events=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_event | cut -d: -f1) ${output}/${Name}.full.tmp)
		kg_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_AF | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		kg_maxfreq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_max_AF | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		#GD_ID=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_ID | cut -d: -f1) ${output}/${Name}.full.tmp)
		#GD_N_HOMALT=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_N_HOMALT | cut -d: -f1) ${output}/${Name}.full.tmp)
		#GD_AT=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_AF | cut -d: -f1) ${output}/${Name}.full.tmp)
		#GD_POPMAX=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_POPMAX_AF | cut -d: -f1) ${output}/${Name}.full.tmp)
		#GD_ID_others=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_ID_others | cut -d: -f1) ${output}/${Name}.full.tmp)
		IMH_ID=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx IMH_ID | cut -d: -f1) ${output}/${Name}.full.tmp)
		IMH_AF=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx IMH_AF | cut -d: -f1) ${output}/${Name}.full.tmp)
		IMH_ID_others=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx IMH_ID_others | cut -d: -f1) ${output}/${Name}.full.tmp)
		dbvar_events=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_event | cut -d: -f1) ${output}/${Name}.full.tmp)
		dbvar_IDs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_variant | cut -d: -f1) ${output}/${Name}.full.tmp)
		dbvar_status=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_status | cut -d: -f1) ${output}/${Name}.full.tmp)
		TADcoordinates=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx TADcoordinates | cut -d: -f1) ${output}/${Name}.full.tmp)
		ENCODE_expiriments=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx ENCODEexperiments | cut -d: -f1) ${output}/${Name}.full.tmp)
		GC_left=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GCcontent_left | cut -d: -f1) ${output}/${Name}.full.tmp)
		GC_right=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GCcontent_right | cut -d: -f1) ${output}/${Name}.full.tmp)
		repeat_loc_left=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Repeats_coord_left | cut -d: -f1) ${output}/${Name}.full.tmp)
		repeat_type_left=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Repeats_type_left | cut -d: -f1) ${output}/${Name}.full.tmp)
		repeat_loc_right=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Repeats_coord_right | cut -d: -f1) ${output}/${Name}.full.tmp)
		repeat_type_right=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Repeats_type_right | cut -d: -f1) ${output}/${Name}.full.tmp)
		gene=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Gene_name | cut -d: -f1) ${output}/${Name}.split.tmp))
		#nm=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx NM | cut -d: -f1) ${output}/${Name}.split.tmp))
		CDS_len=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx CDS_length | cut -d: -f1) ${output}/${Name}.split.tmp))
		tx_len=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx tx_length | cut -d: -f1) ${output}/${Name}.split.tmp))
		loc=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx location | cut -d: -f1) ${output}/${Name}.split.tmp))
		int_st=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx intersectStart | cut -d: -f1) ${output}/${Name}.split.tmp))
		int_end=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx intersectEnd | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_gain_IDs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_IDs | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_gain_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp| sed 's/_//g'))
		sp_DGV_loss_IDs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_IDs | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_loss_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp | sed 's/_//g'))
		#sp_DDD_SVs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_SV | cut -d: -f1) ${output}/${Name}.split.tmp))
		#sp_DDD_gain_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DUP_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp))
		#sp_DDD_loss_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DEL_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_kg_events=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_event | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_kg_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_AF | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_kg_maxfreq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_max_AF | cut -d: -f1) ${output}/${Name}.split.tmp))
		promoter=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx promoters | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_dbvar_events=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_event | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_dbvar_IDs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_variant | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_dbvar_status=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dbVar_status | cut -d: -f1) ${output}/${Name}.split.tmp))
		HI_GC=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx HI_CGscore | cut -d: -f1) ${output}/${Name}.split.tmp))
		TriS_GC=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx TriS_CGscore | cut -d: -f1) ${output}/${Name}.split.tmp))
		DDD_stat=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_status | cut -d: -f1) ${output}/${Name}.split.tmp))
		DDD_mode=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_mode | cut -d: -f1) ${output}/${Name}.split.tmp))
		#DDD_consq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_consequence | cut -d: -f1) ${output}/${Name}.split.tmp)) #NOT REPORTED
		DDD_disease=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_disease | cut -d: -f1) ${output}/${Name}.split.tmp))
		DDD_pmid=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_pmids | cut -d: -f1) ${output}/${Name}.split.tmp))
		DDD_HI=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx HI_DDDpercent | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_synZ=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx synZ_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_misZ=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx misZ_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_pLI=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx pLI_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_delz=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx delZ_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_dupz=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx dupZ_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		ExAC_cnvz=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx cnvZ_ExAC | cut -d: -f1) ${output}/${Name}.split.tmp))
		morbidGene=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx morbidGenes | cut -d: -f1) ${output}/${Name}.split.tmp))
		morbidGene_candidate=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx morbidGenesCandidates | cut -d: -f1) ${output}/${Name}.split.tmp))
		MIM=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Mim_Number | cut -d: -f1) ${output}/${Name}.split.tmp))
		Phenotypes=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Phenotypes | cut -d: -f1) ${output}/${Name}.split.tmp))
		Inheritance=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx Inheritance | cut -d: -f1) ${output}/${Name}.split.tmp))
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${Freq}\t${DGV_gain_IDs}\t${DGV_gain_freq}\t${DGV_loss_IDs}\t${DGV_loss_freq}\t${sp_DGV_gain_IDs}\t${sp_DGV_gain_freq}\t${sp_DGV_loss_IDs}\t${sp_DGV_loss_freq}\t${kg_events}\t${kg_freq}\t${kg_maxfreq}\t${sp_kg_events}\t${sp_kg_freq}\t${sp_kg_maxfreq}\t$IMH_ID\t$IMH_AF\t$IMH_ID_others\t${dbvar_events}\t${dbvar_IDs}\t${dbvar_status}\t${sp_dbvar_events}\t${sp_dbvar_IDs}\t${sp_dbvar_status}\t${TADcoordinates}\t${ENCODE_expiriments}\t${GC_left}\t${GC_right}\t${repeat_loc_left}\t${repeat_type_left}\t${repeat_loc_right}\t${repeat_type_right}\t${gene}\t${CDS_len}\t${tx_len}\t${loc}\t${int_st}\t${int_end}\t${promoter}\t${HI_GC}\t${TriS_GC}\t${DDD_stat}\t${DDD_mode}\t${DDD_disease}\t${DDD_pmid}\t${DDD_HI}\t${ExAC_synZ}\t${ExAC_misZ}\t${ExAC_pLI}\t${ExAC_delz}\t${ExAC_dupz}\t${ExAC_cnvz}\t${morbidGene}\t${morbidGene_candidate}\t${MIM}\t${Phenotypes//,}\t${Inheritance//,}" >> ${output}/ALL.sort.ann.bed
		rm ${output}/${Name}.full.tmp
		rm ${output}/${Name}.split.tmp
	done < ${output}/ALL.sort.unann.bed
fi
rm ${output}/ALL.sort.unann.bed
rm -r ${output}/Annotat 

####################################################
echo -e "STATUS [$(date)]: Performing Custom Annotation..."
####################################################


echo -e "STATUS [$(date)]: Annotating OMIM Symptoms for Genes and Promoters..."
python ${annotation_dir}/Custom_omim_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Genes_to_phenotypes.txt > ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating HPA data..."
python ${annotation_dir}/Custom_hpa_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/HPA > ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating GTEx data..."
python ${annotation_dir}/Custom_gtex_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/GTEx/GTEx_MAX_TISSUE.gct 1 > ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating Brainspan data..."
python ${annotation_dir}/Custom_brainspan_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Brainspan 0.5 > ${output}/ALL.sort.ann.bed2 ##Getting syntax error
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed ### ALL COMING BACK AS NONE?
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating ENCODE TF data..."
bash ${annotation_dir}/Custom_Reg_ann_38.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data > ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2

####################################################
echo -e "STATUS [$(date)]: Hard Filtering..."
####################################################

while read Chr Start Stop Type Name Score Filter InternalF Oth ; do
	if [ ${Type} = "Type" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}\tFail_Reason" > ${output}/filtFAIL.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" > ${output}/All.filt.bed
	elif [ ${Type} = "DEL" ]; then
		DGV_f=$(echo -e "${Oth}" | cut -f4 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		DDD_f=$(echo -e "${Oth}" | cut -f11 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		gN_f=$(echo -e "${Oth}" | cut -f24 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		if [ $(echo -e "${Oth}" | cut -f18) = "DEL" ] || [ $(echo -e "${Oth}" | cut -f15) = "DEL;DUP" ]; then
			kg_f=$(echo -e "${Oth}" | cut -f17 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		else
			kg_f="Pass"
		fi
		if [ ${DDD_f} = "Fail" ] || [ ${DGV_f} = "Fail" ] || [ ${kg_f} = "Fail" ]; then
			f_reason="Frequency_in"
			if [ ${DDD_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_DDD")
			fi
			if [ ${DGV_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_DGV")
			fi
			if [ ${kg_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_dbVar")
			fi
			if [ ${gN_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_gNomad")
			fi
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}\t${f_reason}" >> ${output}/filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/All.filt.bed
		fi
	elif [ ${Type} = "DUP" ] ; then
		DDD_f=$(echo -e "${Oth}" | cut -f10 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		DGV_f=$(echo -e "${Oth}" | cut -f1 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		gN_f=$(echo -e "${Oth}" | cut -f24 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		if [ $(echo -e "${Oth}" | cut -f18) = "DUP" ] || [ $(echo -e "${Oth}" | cut -f15) = "DEL;DUP" ]; then
			kg_f=$(echo -e "${Oth}" | cut -f17 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
		else
			kg_f="Pass"
		fi
		if [ ${DDD_f} = "Fail" ] || [ ${DGV_f} = "Fail" ] || [ ${kg_f} = "Fail" ]; then
			f_reason="Frequency_in"
			if [ ${DDD_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_DDD")
			fi
			if [ ${DGV_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_DGV")
			fi
			if [ ${kg_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_dbVar")
			fi
			if [ ${gN_f} = "Fail" ]; then
				f_reason=$(echo -e "${f_reason}_gNomad")
			fi
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}\t${f_reason}" >> ${output}/filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/All.filt.bed
		fi
	else
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/All.filt.bed
	fi
done < ${output}/ALL.sort.ann.bed
#rm ${output}/ALL.sort.ann.bed

#######################
echo -e "STATUS [$(date)]: Merging Records..."
#######################
if [ -e ${output}/TRLXa.bed ]; then
	rm ${output}/TRLXa.bed
fi
if [ -e ${output}/TRLXb.bed ]; then
	rm ${output}/TRLXb.bed
fi
while read Chr Start Stop Type Name Score Filter InternalF Oth ; do
	if [ ${Type} = 'Type' ]; then
		echo -e "${Chr}_A\t${Start}_A\t${Stop}_A\tChr_B\tStart_B\tStop_B\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" > ${output}/nonTRLX.bedpe
		echo -e "${Chr}_A\t${Start}_A\t${Stop}_A\tChr_B\tStart_B\tStop_B\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" > ${output}/TRLX.bedpe
	elif [ ${Type} = "TRLX" ] || [ ${Type} = 'ORPHAN_TRLX' ]; then
		if [[ ${Name} == *A ]]; then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/TRLXa.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/TRLXb.bed
		fi
	else
		echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" >> ${output}/nonTRLX.bedpe
	fi
done < ${output}/All.filt.bed
while read Chr Start Stop Type Name Score Filter InternalF Oth ; do
	suffix="_A"
	Name=${Name%"$suffix"}
	grep "${Name}" ${output}/TRLXb.bed > ${output}/TRLX_tmpB.bed
	echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Score}\t${Filter}\t${InternalF}\t${Oth}" > ${output}/TRLX_tmpA.bed
	python ${annotation_dir}/Merge_line.py ${output}/TRLX_tmpA.bed ${output}/TRLX_tmpB.bed >> ${output}/TRLX.bedpe
done < ${output}/TRLXa.bed

head -n 1 ${output}/TRLX.bedpe > ${output}/ALL.sort.bedpe
cat ${output}/nonTRLX.bedpe ${output}/TRLX.bedpe > ${output}/ALL.unsort.bedpe

for chr in $( seq 1 22 ) X Y; do
	awk -v OFS="\t" -v CHR=$chr '{if ($1==CHR) print}' ${output}/ALL.unsort.bedpe |sort -k1,1n -k2,2n -k3,3n >> ${output}/ALL.sort.bedpe
done


#######################
echo -e "STATUS [$(date)]: Cleaning Up..."
rm ${output}/ALL.unsort.bedpe ${output}/nonTRLX.bedpe ${output}/TRLX.bedpe 
rm ${output}/TRLX_tmpA.bed ${output}/TRLX_tmpB.bed ${output}/TRLXa.bed  ${output}/TRLXb.bed
rm ${output}/All.filt.bed ${output}/Temp.Match.bed
rm ${output}/TRLX.bed ${output}/INV.bed ${output}/DUP.bed ${output}/DEL.bed ${output}/ORPHAN_TRLX.bed
#######################
