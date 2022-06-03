#!/bin/bash

#################################
# Converter and annotation      #
# script for use with HOLMES    #
#################################

output=$1
params=$2
sample_list=$3
freq=$4
max_samp=$5
annotation_dir=$6

min_samples=1000 #min samples for DGV freq

####################################################
echo -e "STATUS [$(date)]: Starting..."
####################################################

module load BEDTools/2.26.0-gimkl-2017a
module load Tcl/8.6.6-gimkl-2017a
. ${params}
#Source params file

function join_by { local IFS="$1"; shift; echo "$*"; }

if ! [ -e ${output} ]; then
	mkdir ${output}
fi
if ! [ -e ${output}/Annotat ]; then
	mkdir ${output}/Annotat
fi

####################################################
echo -e "STATUS [$(date)]: Making master BED file..."
####################################################

awk -v OFS="\t" '{ print $1, $2, $3, "DUP", $4 }' ${OUTDIR}/SV_calls/*.duplication.bed | fgrep -v "#" > ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $1, $2, $3, "DEL", $4 }' ${OUTDIR}/SV_calls/*.deletion.bed | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
while read Chr Start Stop Name Type Obs Samples Clusters; do
	if [ ${Chr} != "#chr" ]; then
		if [ ${Type} != "deletion" ] && [ ${Type} != "duplication" ]; then ##here
			echo -e "${Clusters}\t${Name}" | sed 's/\[//g' | sed 's/\]//g' >> ${output}/Unres.mult.ref.txt
			echo -e "${Clusters}" | sed 's/\[//g' | sed 's/\]//g' >> ${output}/Unres.mult.list.txt
		else
			echo -e "${Chr}\t${Start}\t${Stop}\tUNR\t${Name}" | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
		fi
	fi
done < ${OUTDIR}/SV_calls/*.unresolved.bed
sed -i 's/,/\n/g' ${output}/Unres.mult.list.txt
while read Cluster; do
	typ=$(echo -e "${Cluster}" | cut -d_ -f2)
	awk -v OFS="\t" -v clust=${Cluster} '{ if ($7==clust) print $1, $2, $3, "UNR", $7"_A" }' ${OUTDIR}/data/clusters/*.${typ}.bkpts.bedpe >> ${output}/ALL.unsort.unann.bed
	awk -v OFS="\t" -v clust=${Cluster} '{ if ($7==clust) print $4, $5, $6, "UNR", $7"_B" }' ${OUTDIR}/data/clusters/*.${typ}.bkpts.bedpe >> ${output}/ALL.unsort.unann.bed
done < ${output}/Unres.mult.list.txt
awk -v OFS="\t" '{ print $1, $2, $3, "TRLX", $7"_A" }' ${OUTDIR}/SV_calls/*.translocation.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $4, $5, $6, "TRLX", $7"_B" }' ${OUTDIR}/SV_calls/*.translocation.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $1, $2, $3, "INS", $7"_A" }' ${OUTDIR}/SV_calls/*.insertion.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $4, $5, $6, "INS", $7"_B" }' ${OUTDIR}/SV_calls/*.insertion.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $1, $2, $3, "INV", $7"_A" }' ${OUTDIR}/SV_calls/*.inversion.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $4, $5, $6, "INV", $7"_B" }' ${OUTDIR}/SV_calls/*.inversion.bedpe | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
awk -v OFS="\t" '{ print $1, $2, $3, $5, $4 }' ${OUTDIR}/SV_calls/*.complex.bed | fgrep -v "#" >> ${output}/ALL.unsort.unann.bed
if [ -e ${output}/ALL.sort.unann.bed ]; then
	rm ${output}/ALL.sort.unann.bed
fi
for chr in $( seq 1 22 ) X Y; do
	awk -v OFS="\t" -v CHR=$chr '{if ($1==CHR) print}' ${output}/ALL.unsort.unann.bed |sort -k1,1n -k2,2n -k3,3n >> ${output}/ALL.sort.unann.bed
done
rm ${output}/ALL.unsort.unann.bed
rm ${output}/Unres.mult.list.txt


####################################################
echo -e "STATUS [$(date)]: Running AnnotSV..."
####################################################

export ANNOTSV=/scale_wlg_persistent/filesets/project/uoa02608/modules/AnnotSV_2.1
$ANNOTSV/bin/AnnotSV -SVinputFile ${output}/ALL.sort.unann.bed -outputdir ${output}/Annotat > ${output}/AnnotatSV.log
sub="$(fgrep -v "#" ${output}/ALL.sort.unann.bed|wc -l|head -n 1)"
ann="$(awk -v OFS="\t" -v type="full" '{if ($8==type) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv|wc -l |head -n 1)"
####################################################
echo -e "STATUS [$(date)]: Compiling annotation..."
####################################################

if [ ${sub} != ${ann} ]; then
	echo -e "ERROR:: annotated records does not equal submitted records: ${sub} , ${ann}"
	exit
else
	echo -e "#Chr\tStart\tStop\tType\tName\tDGV_gain_IDs\tDGV_gain_freq\tDGV_loss_IDs\tDGV_loss_freq\tsp_DGV_gain_IDs\tsp_DGV_gain_freq\tsp_DGV_loss_IDs\tsp_DGV_loss_freq\tDDD_SVs\tDDD_gain_freq\tDDD_loss_freq\tsp_DDD_SVs\tsp_DDD_gain_freq\tsp_DDD_loss_freq\tkg_events\tkg_freq\tkg_maxfreq\tsp_kg_events\tsp_kg_freq\tsp_kg_maxfreq\tGD_ID\tGD_N_HOM\tGD_AT\tGD_POPMAX\tGD_ID_others\tIMH_ID\tIMH_AF\tIMH_ID_others\tdbvar_events\tdbvar_IDs\tdbvar_status\tsp_dbvar_events\tsp_dbvar_IDs\tsp_dbvar_status\tTADcoordinates\tENCODE_expiriments\tGC_left\tGC_right\trepeat_loc_left\trepeat_type_left\trepeat_loc_right\trepeat_type_right\tgene\tnm\tCDS_len\ttx_len\tloc\tint_st\tint_end\tpromoter\tHI_GC\tTriS_GC\tDDD_stat\tDDD_mode\tDDD_disease\tDDD_pmid\tDDD_HI\tExAC_synZ\tExAC_misZ\tExAC_pLI\tExAC_delz\tExAC_dupz\tExAC_cnvz\tmorbidGene\tmorbidGene_candidate\tMIM\tPhenotypes\tInheritance" > ${output}/ALL.sort.ann.bed
	while read Chr Start Stop Type Name ; do
		awk -v OFS="\t" -v type="split" -v chr=${Chr} -v stop=${Stop} -v start=${Start} -v t=${Type} '{if ($6==t && $8==type && $2==chr && $3==start && $4==stop) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/\t\t/\t-\t/g'| sed 's/\t\t/\t-\t/g'| sed 's/ /_/g'| sed 's/\t-1/\t-/g' | sort -u > ${output}/${Name}.split.tmp
		awk -v OFS="\t" -v type="full" -v chr=${Chr} -v stop=${Stop} -v start=${Start} -v t=${Type} '{if ($6==t && $8==type && $2==chr && $3==start && $4==stop) print $0}' ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/\t\t/\t-\t/g'| sed 's/\t\t/\t-\t/g'| sed 's/ /_/g'| sed 's/\t-1/\t-/g' | sort -u > ${output}/${Name}.full.tmp ###HEAD HERE IF ISSUE WITH ANNOTAT CONTINUES
		if [ $(wc -l < ${output}/${Name}.split.tmp | head -n 1) == "0" ]; then
			echo -e "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-" > ${output}/${Name}.split.tmp
		fi
		DGV_gain_IDs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_IDs | cut -d: -f1) ${output}/${Name}.full.tmp)
		DGV_gain_freq=$( cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_n_samples_tested | cut -d: -f1),$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | awk -v min=${min_samples} '{if ( $1>=min ) print $2 ; else print "Fewer_than_"min"_samples"}' | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		DGV_loss_IDs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_IDs | cut -d: -f1) ${output}/${Name}.full.tmp)
		DGV_loss_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_n_samples_tested | cut -d: -f1),$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | awk -v OFS="\t"  -v min=${min_samples} '{if ( $1>=min ) print $2 ; else print "Fewer_than_"min"_samples"}' | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		DDD_SVs=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_SV | cut -d: -f1) ${output}/${Name}.full.tmp)
		DDD_gain_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DUP_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		DDD_loss_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DEL_Frequency | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		kg_events=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_event | cut -d: -f1) ${output}/${Name}.full.tmp)
		kg_freq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_AF | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		kg_maxfreq=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx 1000g_max_AF | cut -d: -f1) ${output}/${Name}.full.tmp | sed 's/_//g' | sed 's/Fewerthan/Fewer_than_/g' | sed 's/samples/_samples/g')
		GD_ID=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_ID | cut -d: -f1) ${output}/${Name}.full.tmp)
		GD_N_HOMALT=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_N_HOMALT | cut -d: -f1) ${output}/${Name}.full.tmp)
		GD_AT=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_AF | cut -d: -f1) ${output}/${Name}.full.tmp)
		GD_POPMAX=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_POPMAX_AF | cut -d: -f1) ${output}/${Name}.full.tmp)
		GD_ID_others=$(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx GD_ID_others | cut -d: -f1) ${output}/${Name}.full.tmp)
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
		nm=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx NM | cut -d: -f1) ${output}/${Name}.split.tmp))
		CDS_len=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx CDS_length | cut -d: -f1) ${output}/${Name}.split.tmp))
		tx_len=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx tx_length | cut -d: -f1) ${output}/${Name}.split.tmp))
		loc=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx location | cut -d: -f1) ${output}/${Name}.split.tmp))
		int_st=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx intersectStart | cut -d: -f1) ${output}/${Name}.split.tmp))
		int_end=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx intersectEnd | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_gain_IDs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_IDs | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_gain_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_GAIN_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp| sed 's/_//g'))
		sp_DGV_loss_IDs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_IDs | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DGV_loss_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DGV_LOSS_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp | sed 's/_//g'))
		sp_DDD_SVs=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_SV | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DDD_gain_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DUP_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp))
		sp_DDD_loss_freq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_DEL_Frequency | cut -d: -f1) ${output}/${Name}.split.tmp))
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
		DDD_consq=$( join_by : $(cut -f$(head -n 1 ${output}/Annotat/ALL.sort.unann.annotated.tsv | sed 's/ /_/g' | sed 's/\t/\n/g' | grep -nx DDD_consequence | cut -d: -f1) ${output}/${Name}.split.tmp))
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
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${DGV_gain_IDs}\t${DGV_gain_freq}\t${DGV_loss_IDs}\t${DGV_loss_freq}\t${sp_DGV_gain_IDs}\t${sp_DGV_gain_freq}\t${sp_DGV_loss_IDs}\t${sp_DGV_loss_freq}\t${DDD_SVs}\t${DDD_gain_freq}\t${DDD_loss_freq}\t${sp_DDD_SVs}\t${sp_DDD_gain_freq}\t${sp_DDD_loss_freq}\t${kg_events}\t${kg_freq}\t${kg_maxfreq}\t${sp_kg_events}\t${sp_kg_freq}\t${sp_kg_maxfreq}\t$GD_ID\t$GD_N_HOMALT\t$GD_AT\t$GD_POPMAX\t$GD_ID_others\t$IMH_ID\t$IMH_AF\t$IMH_ID_others\t${dbvar_events}\t${dbvar_IDs}\t${dbvar_status}\t${sp_dbvar_events}\t${sp_dbvar_IDs}\t${sp_dbvar_status}\t${TADcoordinates}\t${ENCODE_expiriments}\t${GC_left}\t${GC_right}\t${repeat_loc_left}\t${repeat_type_left}\t${repeat_loc_right}\t${repeat_type_right}\t${gene}\t${nm}\t${CDS_len}\t${tx_len}\t${loc}\t${int_st}\t${int_end}\t${promoter}\t${HI_GC}\t${TriS_GC}\t${DDD_stat}\t${DDD_mode}\t${DDD_disease}\t${DDD_pmid}\t${DDD_HI}\t${ExAC_synZ}\t${ExAC_misZ}\t${ExAC_pLI}\t${ExAC_delz}\t${ExAC_dupz}\t${ExAC_cnvz}\t${morbidGene}\t${morbidGene_candidate}\t${MIM}\t${Phenotypes}\t${Inheritance}" >> ${output}/ALL.sort.ann.bed
		rm ${output}/${Name}.full.tmp
		rm ${output}/${Name}.split.tmp
	done < ${output}/ALL.sort.unann.bed
fi
rm ${output}/ALL.sort.unann.bed
rm -r ${output}/Annotat 
####################################################
echo -e "STATUS [$(date)]: Performing Custom Annotation..."
####################################################
module load Python/2.7.14-gimkl-2017a
echo -e "STATUS [$(date)]: Performing Custom Annotation..."
echo -e "STATUS [$(date)]: Annotating OMIM Symptoms for Genes and Promoters..."
#bash ${annotation_dir}/Custom_omim_ann.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Genes_to_phenotypes.txt > ${output}/ALL.sort.ann.bed2
python ${annotation_dir}/Custom_omim_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Genes_to_phenotypes.txt > ${output}/ALL.sort.ann.bed2
#head -n 2 ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating HPA data..."
#bash ${annotation_dir}/Custom_hpa_ann.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/HPA > ${output}/ALL.sort.ann.bed2
python ${annotation_dir}/Custom_hpa_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/HPA > ${output}/ALL.sort.ann.bed2
#head -n 2 ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating GTEx data..."
#bash ${annotation_dir}/Custom_gtes_ann.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/GTEx/GTEx_MAX_TISSUE.gct 1 > ${output}/ALL.sort.ann.bed2
python ${annotation_dir}/Custom_gtex_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/GTEx/GTEx_MAX_TISSUE.gct 1 > ${output}/ALL.sort.ann.bed2
#head -n 2 ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating Brainspan data..."
#bash ${annotation_dir}/Custom_brainspan_ann.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Brainspan 0.5 > ${output}/ALL.sort.ann.bed2
python ${annotation_dir}/Custom_brainspan_ann.py ${output}/ALL.sort.ann.bed ${annotation_dir}/Data/Brainspan 0.5 > ${output}/ALL.sort.ann.bed2
#head -n 2 ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed2
echo -e "STATUS [$(date)]: Annotating ENCODE TF data..."
bash ${annotation_dir}/Custom_Reg_ann.sh ${output}/ALL.sort.ann.bed ${annotation_dir}/Data > ${output}/ALL.sort.ann.bed2
#head -n 2 ${output}/ALL.sort.ann.bed2
cp ${output}/ALL.sort.ann.bed2 ${output}/ALL.sort.ann.bed
#rm ${output}/ALL.sort.ann.bed2
####################################################
echo -e "STATUS [$(date)]: Hard Filtering..."
####################################################

while read Chr Start Stop Type Name Oth ; do
	if [ ${Type} = "Type" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/DEL.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/DUP.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/INS.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/INV.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/TRLX.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/UNR.prefilt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/COMP.prefilt.bed
	elif [ ${Type} = "DEL" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/DEL.prefilt.bed
	elif [ ${Type} = "DUP" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/DUP.prefilt.bed
	elif [ ${Type} = "INS" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/INS.prefilt.bed
	elif [ ${Type} = "INV" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/INV.prefilt.bed
	elif [ ${Type} = "TRLX" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/TRLX.prefilt.bed
	elif [ ${Type} = "UNR" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/UNR.prefilt.bed
	else
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/COMP.prefilt.bed
	fi
done < ${output}/ALL.sort.ann.bed
rm ${output}/ALL.sort.ann.bed
echo -e "STATUS [$(date)]: --For Deletions..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/DEL.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/DEL.filtFAIL.bed
	else
		Con=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $6}' ${OUTDIR}/SV_calls/*.deletion.bed)
		if [ ${Con} = "LOW" ]; then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tLOW_confidence_call" >> ${output}/DEL.filtFAIL.bed
		else
			DDD_f=$(echo -e "${Oth}" | cut -f11 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			DGV_f=$(echo -e "${Oth}" | cut -f4 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			gN_f=$(echo -e "${Oth}" | cut -f24 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			if [ $(echo -e "${Oth}" | cut -f15) = "DEL" ] || [ $(echo -e "${Oth}" | cut -f15) = "DEL;DUP" ]; then
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
				echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\t${f_reason}" >> ${output}/DEL.filtFAIL.bed
			else
				Obs=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $7}' ${OUTDIR}/SV_calls/*.deletion.bed | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
				if [ ${Obs} = "Fail" ];then
					echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/DEL.filtFAIL.bed
				else
					echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/DEL.filt.bed
				fi
			fi
		fi
	fi
done < ${output}/DEL.prefilt.bed
rm ${output}/DEL.prefilt.bed
echo -e "STATUS [$(date)]: --For Duplications..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/DUP.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/DUP.filtFAIL.bed
	else
		Con=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $6}' ${OUTDIR}/SV_calls/*.duplication.bed)
		if [ ${Con} = "LOW" ]; then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tLOW_confidence_call" >> ${output}/DUP.filtFAIL.bed
		else
			DDD_f=$(echo -e "${Oth}" | cut -f10 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			DGV_f=$(echo -e "${Oth}" | cut -f1 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			gN_f=$(echo -e "${Oth}" | cut -f24 | awk -v f=${freq} '{if ( $1>=f && $1 ~ /^[0-9/.]+$/ ) print "Fail"; else print "Pass"}')
			if [ $(echo -e "${Oth}" | cut -f15) = "DUP" ] || [ $(echo -e "${Oth}" | cut -f15) = "DEL;DUP" ]; then
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
				echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\t${f_reason}" >> ${output}/DUP.filtFAIL.bed
			else
				Obs=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $7}' ${OUTDIR}/SV_calls/*.duplication.bed | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
				if [ ${Obs} = "Fail" ];then
					echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/DUP.filtFAIL.bed
				else
					echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/DUP.filt.bed
				fi
			fi
		fi
	fi
done < ${output}/DUP.prefilt.bed
rm ${output}/DUP.prefilt.bed
echo -e "STATUS [$(date)]: --For Insertions..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/INS.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/INS.filtFAIL.bed
	else
		Obs=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $10}' ${OUTDIR}/SV_calls/*.insertion.bedpe | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
		if [ ${Obs} = "Fail" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/INS.filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/INS.filt.bed
		fi
	fi
done < ${output}/INS.prefilt.bed
rm ${output}/INS.prefilt.bed
echo -e "STATUS [$(date)]: --For Inversions..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/INV.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/INV.filtFAIL.bed
	else
		Obs=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $10}' ${OUTDIR}/SV_calls/*.inversion.bedpe | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
		if [ ${Obs} = "Fail" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/INV.filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/INV.filt.bed
		fi
	fi
done < ${output}/INV.prefilt.bed
rm ${output}/INV.prefilt.bed
echo -e "STATUS [$(date)]: --For Translocation..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/TRLX.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/TRLX.filtFAIL.bed
	else
		Obs=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $10}' ${OUTDIR}/SV_calls/*.translocation.bedpe | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
		if [ ${Obs} = "Fail" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/TRLX.filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/TRLX.filt.bed
		fi
	fi
done < ${output}/TRLX.prefilt.bed
rm ${output}/TRLX.prefilt.bed
echo -e "STATUS [$(date)]: --For Complex..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/COMP.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/COMP.filtFAIL.bed
	else
		Obs=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $6}' ${OUTDIR}/SV_calls/*.complex.bed | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
		if [ ${Obs} = "Fail" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/COMP.filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/COMP.filt.bed
		fi
	fi
done < ${output}/COMP.prefilt.bed
rm ${output}/COMP.prefilt.bed
echo -e "STATUS [$(date)]: --For Unresolved..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} = "#Chr" ]; then
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/UNR.filt.bed
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}\tFail_Cause" > ${output}/UNR.filtFAIL.bed
	else
		N2=${Name}
		if [ $(echo -e "${Name}" | cut -d_ -f2) != "unresolved" ];then
			Name=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($1==name) print $2}' ${output}/Unres.mult.ref.txt)
		fi
		Obs=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $6}' ${OUTDIR}/SV_calls/*.unresolved.bed | awk -v max=${max_samp} '{if ($1>=max) print "Fail"; else print "Pass"}')
		if [ ${Obs} = "Fail" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${N2}\t${Oth}\tIn_${max_samp}_or_more_samples" >> ${output}/UNR.filtFAIL.bed
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${N2}\t${Oth}" >> ${output}/UNR.filt.bed
		fi
	fi
done < ${output}/UNR.prefilt.bed
rm ${output}/UNR.prefilt.bed
echo -e "STATUS [$(date)]: HARD Breakpoint Filtering Complete...."
echo -e "Deletions Passing: $(fgrep -v "#" ${output}/DEL.filt.bed|wc -l|head -n 1), Deletions Failing: $(fgrep -v "#" ${output}/DEL.filtFAIL.bed|wc -l|head -n 1)" > ${output}/HARDfiltering.log
echo -e "Duplications Passing: $(fgrep -v "#" ${output}/DUP.filt.bed|wc -l|head -n 1), Duplications Failing: $(fgrep -v "#" ${output}/DUP.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
echo -e "Insertions Passing: $(fgrep -v "#" ${output}/INS.filt.bed|wc -l|head -n 1), Insertions Failing: $(fgrep -v "#" ${output}/INS.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
echo -e "Inversions Passing: $(fgrep -v "#" ${output}/INV.filt.bed|wc -l|head -n 1), Inversions Failing: $(fgrep -v "#" ${output}/INV.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
echo -e "Translocations Passing: $(fgrep -v "#" ${output}/TRLX.filt.bed|wc -l|head -n 1), Translocations Failing: $(fgrep -v "#" ${output}/TRLX.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
echo -e "Complex Passing: $(fgrep -v "#" ${output}/COMP.filt.bed|wc -l|head -n 1), Complex Failing: $(fgrep -v "#" ${output}/COMP.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
echo -e "Unresolved Passing: $(fgrep -v "#" ${output}/UNR.filt.bed|wc -l|head -n 1), Unresolved Failing: $(fgrep -v "#" ${output}/UNR.filtFAIL.bed|wc -l|head -n 1)" >> ${output}/HARDfiltering.log
#######################
#Split and report
echo -e "STATUS [$(date)]: Splitting and Reporting By Sample..."
#######################

while read Sample Oth ; do
	echo -e "#Chr_A\tStart_A\tStop_A\tChr_B\tStart_B\tStop_B\tOri\tType\tName\t$(head -n 1 ${output}/DEL.filt.bed | cut -f6-)" > ${output}/${Sample}.vars.bedpe
done < ${sample_list}
echo -e "#Chr_A\tStart_A\tStop_A\tChr_B\tStart_B\tStop_B\tOri\tType\tName\t$(head -n 1 ${output}/DEL.filt.bed | cut -f6-)" > ${output}/OTH.vars.bedpe
echo -e "STATUS [$(date)]: --For Deletions..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} != "#Chr" ]; then
		samps=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $8}' ${OUTDIR}/SV_calls/*.deletion.bed | sed 's/\]//g'| sed 's/\[//g')
		for sam in $(echo $samps | tr "," "\n") ; do
			if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/OTH.vars.bedpe
			else
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/${sam}.vars.bedpe
			fi
		done
	fi
done < ${output}/DEL.filt.bed
echo -e "STATUS [$(date)]: --For Duplications..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} != "#Chr" ]; then
		samps=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $8}' ${OUTDIR}/SV_calls/*.duplication.bed | sed 's/\]//g'| sed 's/\[//g')
		for sam in $(echo $samps | tr "," "\n") ; do
			if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/OTH.vars.bedpe
			else
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/${sam}.vars.bedpe
			fi
		done
	fi
done < ${output}/DUP.filt.bed
echo -e "STATUS [$(date)]: --For Insertions..."
for nam in $(awk '{if ($1!="#Chr")print $5}' ${output}/INS.filt.bed |sed 's/_A//g'| sed 's/_B//g'| sort -u) ; do
	oth=""
	for col in $( seq 6 $(head -n 1 ${output}/INS.filt.bed | awk '{print NF}') ) ; do
		oth="${oth}\t$(grep "${nam}_A" ${output}/INS.filt.bed | cut -f${col})|$(grep "${nam}_B" ${output}/INS.filt.bed | cut -f${col})"
	done
	samps=$(awk -v OFS="\t" -v name=$(echo -e "${nam}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $11}' ${OUTDIR}/SV_calls/*.insertion.bedpe | sed 's/\]//g'| sed 's/\[//g')
	for sam in $(echo $samps | tr "," "\n") ; do
		if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.insertion.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.insertion.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.insertion.bedpe)\tINS\t${nam}${oth}" >> ${output}/OTH.vars.bedpe
		else
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.insertion.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.insertion.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.insertion.bedpe)\tINS\t${nam}${oth}" >> ${output}/${sam}.vars.bedpe
		fi
	done
done
echo -e "STATUS [$(date)]: --For Inversions..."
for nam in $(awk '{if ($1!="#Chr")print $5}' ${output}/INV.filt.bed |sed 's/_A//g'| sed 's/_B//g'| sort -u) ; do
	oth=""
	for col in $( seq 6 $(head -n 1 ${output}/INV.filt.bed | awk '{print NF}') ) ; do
		oth="${oth}\t$(grep "${nam}_A" ${output}/INV.filt.bed | cut -f${col})|$(grep "${nam}_B" ${output}/INV.filt.bed | cut -f${col})"
	done
	samps=$(awk -v OFS="\t" -v name=$(echo -e "${nam}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $11}' ${OUTDIR}/SV_calls/*.inversion.bedpe | sed 's/\]//g'| sed 's/\[//g')
	for sam in $(echo $samps | tr "," "\n") ; do
		if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.inversion.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.inversion.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.inversion.bedpe)\tINV\t${nam}${oth}" >> ${output}/OTH.vars.bedpe
		else
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.inversion.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.inversion.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.inversion.bedpe)\tINV\t${nam}${oth}" >> ${output}/${sam}.vars.bedpe
		fi
	done
done
echo -e "STATUS [$(date)]: --For Translocations..."
for nam in $(awk '{if ($1!="#Chr")print $5}' ${output}/TRLX.filt.bed |sed 's/_A//g'| sed 's/_B//g'| sort -u) ; do
	oth=""
	for col in $( seq 6 $(head -n 1 ${output}/TRLX.filt.bed | awk '{print NF}') ) ; do
		oth="${oth}\t$(grep "${nam}_A" ${output}/TRLX.filt.bed | cut -f${col})|$(grep "${nam}_B" ${output}/TRLX.filt.bed | cut -f${col})"
	done
	samps=$(awk -v OFS="\t" -v name=$(echo -e "${nam}" | sed 's/_A//g'| sed 's/_B//g') '{if ($7==name) print $11}' ${OUTDIR}/SV_calls/*.translocation.bedpe | sed 's/\]//g'| sed 's/\[//g')
	for sam in $(echo $samps | tr "," "\n") ; do
		if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.translocation.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.translocation.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.translocation.bedpe)\tTRLX\t${nam}${oth}" >> ${output}/OTH.vars.bedpe
		else
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/SV_calls/*.translocation.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.translocation.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $8}' ${OUTDIR}/SV_calls/*.translocation.bedpe)\tTRLX\t${nam}${oth}" >> ${output}/${sam}.vars.bedpe
		fi
	done
done
echo -e "STATUS [$(date)]: --For Complex..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} != "#Chr" ]; then
		samps=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $7}' ${OUTDIR}/SV_calls/*.complex.bed | sed 's/\]//g'| sed 's/\[//g')
		for sam in $(echo $samps | tr "," "\n") ; do
			if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/OTH.vars.bedpe
			else
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}\t${Name}\t${Oth}" >> ${output}/${sam}.vars.bedpe
			fi
		done
	fi
done < ${output}/COMP.filt.bed
echo -e "STATUS [$(date)]: --For Unresolved..."
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} != "#Chr" ]; then
		if [ $(echo -e "${Name}" | cut -d_ -f2) != "unresolved" ];then
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/UNR.comp.tmp
		else
			echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" >> ${output}/UNR.simp.tmp
		fi
	else
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/UNR.simp.tmp
		echo -e "${Chr}\t${Start}\t${Stop}\t${Type}\t${Name}\t${Oth}" > ${output}/UNR.comp.tmp
	fi
done < ${output}/UNR.filt.bed
while read Chr Start Stop Type Name Oth ; do
	if [ ${Chr} != "#Chr" ]; then
		samps=$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $7}' ${OUTDIR}/SV_calls/*.unresolved.bed | sed 's/\]//g'| sed 's/\[//g')
		for sam in $(echo $samps | tr "," "\n") ; do
			if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}_$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $5}' ${OUTDIR}/SV_calls/*.unresolved.bed)\t${Name}\t${Oth}" >> ${output}/OTH.vars.bedpe
			else
				echo -e "${Chr}\t${Start}\t${Stop}\t-\t-\t-\t-\t${Type}_$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $5}' ${OUTDIR}/SV_calls/*.unresolved.bed)\t${Name}\t${Oth}" >> ${output}/${sam}.vars.bedpe
			fi
		done
	fi
done < ${output}/UNR.simp.tmp
rm ${output}/UNR.simp.tmp
for nam in $(awk '{if ($1!="#Chr")print $5}' ${output}/UNR.comp.tmp |sed 's/_A//g'| sed 's/_B//g'| sort -u) ; do
	oth=""
	for col in $( seq 6 $(head -n 1 ${output}/UNR.comp.tmp | awk '{print NF}') ) ; do
		oth="${oth}\t$(grep "${nam}_A" ${output}/UNR.comp.tmp | cut -f${col})|$(grep "${nam}_B" ${output}/UNR.comp.tmp | cut -f${col})"
	done
	Name=$(awk -v OFS="\t" -v name=$(echo -e "${nam}" | sed 's/_A//g'| sed 's/_B//g') '{if ($1==name) print $2}' ${output}/Unres.mult.ref.txt)
	samps=$(awk -v OFS="\t" -v name=${Name} '{if ($4==name) print $7}' ${OUTDIR}/SV_calls/*.unresolved.bed | sed 's/\]//g'| sed 's/\[//g')
	for sam in $(echo $samps | tr "," "\n") ; do
		if [ $(echo ">$(echo ${sam} | sed 's/[0-9]//g')") = ">" ] ; then
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $9}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $10}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)\tUNK_$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $5}' ${OUTDIR}/SV_calls/*.unresolved.bed)\t${Name}[${nam}]${oth}" >> ${output}/OTH.vars.bedpe
		else
			echo -e "$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)\t$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $9}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)|$(awk -v OFS="\t" -v name=${nam} '{if ($7==name) print $10}' ${OUTDIR}/data/clusters/*.$(echo -e "${nam}" | cut -d_ -f2).bkpts.bedpe)\tUNK_$(awk -v OFS="\t" -v name=$(echo -e "${Name}" | sed 's/_A//g'| sed 's/_B//g') '{if ($4==name) print $5}' ${OUTDIR}/SV_calls/*.unresolved.bed)\t${Name}[${nam}]${oth}" >> ${output}/${sam}.vars.bedpe
		fi
	done
done
rm ${output}/UNR.comp.tmp
rm ${output}/Unres.mult.ref.txt
for type in DEL DUP INS INV TRLX COMP UNR; do
	rm ${output}/${type}.filt.bed
done
#######################
echo -e "STATUS [$(date)]: Cleaning Up..."
#######################
head -n 1 ${output}/DEL.filtFAIL.bed > ${output}/ALL.filtFAIL.bed
for type in DEL DUP INS INV TRLX COMP UNR; do
	cat ${output}/ALL.filtFAIL.bed ${output}/${type}.filtFAIL.bed > ${output}/ALL.filtFAIL.bed2
	cp ${output}/ALL.filtFAIL.bed2 ${output}/ALL.filtFAIL.bed
	rm ${output}/ALL.filtFAIL.bed2 ${output}/${type}.filtFAIL.bed
done
#rm -r ${output}/
