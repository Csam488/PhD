#!/bin/bash
ann_folder=$2
function join_by { local IFS="$1"; shift; echo "$*"; }

while read line ; do
	if [[ $line =~ ^# ]]; then
		g_field=$(echo $line | sed 's/ /\n/g' | grep -nx gene | cut -d: -f1)
		p_field=$(echo $line | sed 's/ /\n/g' | grep -nx promoter | cut -d: -f1)
		echo -e "${line}\tGene_subcellular_location\tGene_RNA_expression_catigory\tGene_RNA_enriched_tissue(s)\tGene_HIGH_expression_protein_tissue\tGene_MED_protein_expression_tissue\tPromoter_subcellular_location\tPromoter_RNA_expression_catigory\tPromoter_RNA_enriched_tissue(s)\tPromorer_HIGH_expression_protein_tissue\tPromoter_MED_protein_expression_tissue"
	else
		g_sub_loc=''
		g_RNA_cat=''
		g_RNA_tissue=''
		g_high_pro=''
		g_med_pro=''
		p_sub_loc=''
		p_RNA_cat=''
		p_RNA_tissue=''
		p_high_pro=''
		p_med_pro=''
		genes=$(echo $line | cut -d' ' -f$g_field | sed 's/:/\n/g')
		promoters=$(echo $line | cut -d' ' -f$p_field | sed 's/:/\n/g')
		for gene in $genes; do
			if [ $gene != '-' ]; then
				PA=$(fgrep -v '#' $2/proteinatlas.tdt | awk -v gen=$gene '{if ($1 == gen || $2 ~ gen ) print $0}')
				sub_loc=$( join_by / $(echo $PA | cut -d' ' -f13))
				RNA_cat=$( join_by / $(echo $PA | cut -d' ' -f16))
				RNA_tissue=$( join_by / $(echo $PA | cut -d' ' -f18 | sed 's/:/-/g'))
				high_pro=$( join_by / $(fgrep -v '#' $2/normal_tissue.tdt | awk -v gen=$gene '{if ($2 == gen ) print $0}' | awk -v lev='High' '{if ($5 == lev ) print $3"-"$4}'))
				med_pro=$( join_by / $(fgrep -v '#' $2/normal_tissue.tdt | awk -v gen=$gene '{if ($2 == gen ) print $0}' | awk -v lev='Medium' '{if ($5 == lev ) print $3"-"$4}'))
				if [ -z $sub_loc ]; then
					sub_loc='-'
				fi
				if [ -z $RNA_cat ]; then
					RNA_cat='-'
				fi
				if [ -z $RNA_tissue ]; then
					RNA_tissue='-'
				fi
				if [ -z $high_pro ]; then
					high_pro='-'
				fi
				if [ -z $med_pro ]; then
					med_pro='-'
				fi
				g_sub_loc=$(echo -e "${g_sub_loc}:${sub_loc}")
				g_RNA_cat=$(echo -e "${g_RNA_cat}:${RNA_cat}")
				g_RNA_tissue=$(echo -e "${g_RNA_tissue}:${RNA_tissue}")
				g_high_pro=$(echo -e "${g_high_pro}:${high_pro}")
				g_med_pro=$(echo -e "${g_med_pro}:${med_pro}")
			else
				g_sub_loc=$(echo -e "${g_sub_loc}:-")
				g_RNA_cat=$(echo -e "${g_RNA_cat}:-")
				g_RNA_tissue=$(echo -e "${g_RNA_tissue}:-")
				g_high_pro=$(echo -e "${g_high_pro}:-")
				g_med_pro=$(echo -e "${g_med_pro}:-")
			fi
		done
		for promoter in $promoters; do
			if [ $promoter != '-' ]; then
				PA=$(fgrep -v '#' $2/proteinatlas.tdt | awk -v gen=$promoter '{if ($1 == gen || $2 ~ gen ) print $0}')
				sub_loc=$( join_by / $(echo $PA | cut -d' ' -f13))
				RNA_cat=$( join_by / $(echo $PA | cut -d' ' -f16))
				RNA_tissue=$( join_by / $(echo $PA | cut -d' ' -f18 | sed 's/:/-/g'))
				high_pro=$( join_by / $(fgrep -v '#' $2/normal_tissue.tdt | awk -v gen=$promoter '{if ($2 == gen ) print $0}' | awk -v lev='High' '{if ($5 == lev ) print $3"-"$4}'))
				med_pro=$( join_by / $(fgrep -v '#' $2/normal_tissue.tdt | awk -v gen=$promoter '{if ($2 == gen ) print $0}' | awk -v lev='Medium' '{if ($5 == lev ) print $3"-"$4}'))
				if [ -z $sub_loc ]; then
					sub_loc='-'
				fi
				if [ -z $RNA_cat ]; then
					RNA_cat='-'
				fi
				if [ -z $RNA_tissue ]; then
					RNA_tissue='-'
				fi
				if [ -z $high_pro ]; then
					high_pro='-'
				fi
				if [ -z $med_pro ]; then
					med_pro='-'
				fi
				p_sub_loc=$(echo -e "${g_sub_loc}:${sub_loc}")
				p_RNA_cat=$(echo -e "${g_RNA_cat}:${RNA_cat}")
				p_RNA_tissue=$(echo -e "${g_RNA_tissue}:${RNA_tissue}")
				p_high_pro=$(echo -e "${g_high_pro}:${high_pro}")
				p_med_pro=$(echo -e "${g_med_pro}:${med_pro}")
			else
				p_sub_loc=$(echo -e "${g_sub_loc}:-")
				p_RNA_cat=$(echo -e "${g_RNA_cat}:-")
				p_RNA_tissue=$(echo -e "${g_RNA_tissue}:-")
				p_high_pro=$(echo -e "${g_high_pro}:-")
				p_med_pro=$(echo -e "${g_med_pro}:-")
			fi
		done
		g_sub_loc=${g_sub_loc#?}
		g_RNA_cat=${g_RNA_cat#?}
		g_RNA_tissue=${g_RNA_tissue#?}
		g_high_pro=${g_high_pro#?}
		g_med_pro=${g_med_pro#?}
		p_sub_loc=${p_sub_loc#?}
		p_RNA_cat=${p_RNA_cat#?}
		p_RNA_tissue=${p_RNA_tissue#?}
		p_high_pro=${p_high_pro#?}
		p_med_pro=${p_med_pro#?}
		echo -e "${line}\t${g_sub_loc}\t${g_RNA_cat}\t${g_RNA_tissue}\t${g_high_pro}\t${g_med_pro}\t${p_sub_loc}\t${p_RNA_cat}\t${p_RNA_tissue}\t${p_high_pro}\t${p_med_pro}"
	fi
done < $1
