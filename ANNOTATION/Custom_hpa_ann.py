import sys
from sys import argv
import re
script, input, source = argv

genes = {}

with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split()
			gs = col[gene].split(':')
			for g in gs:
				if g != '-':
					genes[g] = 1
			gs = col[promoter].split(':')
			for g in gs:
				if g != '-':
					genes[g] = 1
		else:
			col = line.split()
			gene = col.index('gene')
			promoter = col.index('promoter')

High_Protein = {}
Medium_Protein = {}
RNA_tissue = {}
RNA_category = {}
RNA_subcellular_loc = {}
for g in genes.keys():
	High_Protein[g] = []
	Medium_Protein[g] = []
	RNA_tissue[g] = []
	RNA_category[g] = []
	RNA_subcellular_loc[g] = []

with open(source+"/normal_tissue.tdt", "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split('\t')
			if col[1] in genes.keys():
				if col[4] == 'Medium':
					Medium_Protein[col[1]].append(col[2]+"("+col[3]+")")
				if col[4] == 'High':
					High_Protein[col[1]].append(col[2]+"("+col[3]+")")

with open(source+"/proteinatlas.tdt", "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split('\t')
			if col[0] in genes.keys():
				if col[12] != '-':
					RNA_subcellular_loc[col[0]].append(col[12])
				if col[15] != '-':
					RNA_category[col[0]].append(col[15])
				if col[17] != '-':
					RNA_tissue[col[0]].append(col[17])

for g in genes.keys():
	if len(High_Protein[g]) == 0:
		High_Protein[g] = '-'
	else:
		High_Protein[g] = ','.join(High_Protein[g])
	if len(Medium_Protein[g]) == 0:
		Medium_Protein[g] = '-'
	else:
		Medium_Protein[g] = ','.join(Medium_Protein[g])
	if len(RNA_tissue[g]) == 0:
		RNA_tissue[g] = '-'
	else:
		RNA_tissue[g] = ','.join(RNA_tissue[g])
	if len(RNA_category[g]) == 0:
		RNA_category[g] = '-'
	else:
		RNA_category[g] = ','.join(RNA_category[g])
	if len(RNA_subcellular_loc[g]) == 0:
		RNA_subcellular_loc[g] = '-'
	else:
		RNA_subcellular_loc[g] = ','.join(RNA_subcellular_loc[g])
High_Protein['-'] = '-'
Medium_Protein['-'] = '-'
RNA_category['-'] = '-'
RNA_subcellular_loc['-'] = '-'
RNA_tissue['-'] = '-'


with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split()
			gs = col[gene].split(':')
			g_HIGH = []
			g_MED = []
			g_RNAC = []
			g_RNAS = []
			g_RNAT = []
			for g in gs:
				g_HIGH.append(High_Protein[g])
				g_MED.append(Medium_Protein[g])
				g_RNAC.append(RNA_category[g])
				g_RNAS.append(RNA_subcellular_loc[g])
				g_RNAT.append(RNA_tissue[g])
			ps = col[promoter].split(':')
			p_HIGH = []
			p_MED = []
			p_RNAC = []
			p_RNAS = []
			p_RNAT = []
			for p in ps:
				p_HIGH.append(High_Protein[p])
				p_MED.append(Medium_Protein[p])
				p_RNAC.append(RNA_category[p])
				p_RNAS.append(RNA_subcellular_loc[p])
				p_RNAT.append(RNA_tissue[p])
			print line.strip('\n')+'\t'+':'.join(g_RNAS)+'\t'+':'.join(g_RNAC)+'\t'+':'.join(g_RNAT)+'\t'+':'.join(g_HIGH)+'\t'+':'.join(g_MED)+'\t'+':'.join(p_RNAS)+'\t'+':'.join(p_RNAC)+'\t'+':'.join(p_RNAT)+'\t'+':'.join(p_HIGH)+'\t'+':'.join(p_MED)
		else:
			print line.strip('\n')+"\tGene_subcellular_location\tGene_RNA_expression_catigory\tGene_RNA_enriched_tissue(s)\tGene_HIGH_expression_protein_tissue\tGene_MED_protein_expression_tissue\tPromoter_subcellular_location\tPromoter_RNA_expression_catigory\tPromoter_RNA_enriched_tissue(s)\tPromorer_HIGH_expression_protein_tissue\tPromoter_MED_protein_expression_tissue"
'''
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
'''