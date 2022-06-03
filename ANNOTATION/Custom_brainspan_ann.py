import sys
from sys import argv
import re
script, input, source, thresh = argv

genes = {}

with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.strip('\n').split('\t')
			gs = col[gene].split(':')
			for g in gs:
				if g != '-':
					genes[g] = 1
			gs = col[promoter].split(':')
			for g in gs:
				if g != '-':
					genes[g] = 1
		else:
			col = line.strip('\n').split('\t')
			gene = col.index('gene')
			promoter = col.index('promoter')

period = {}
structure = {}
for g in genes.keys():
	period[g] = []
	structure[g] = []

with open(source+'/Period_ref.tdt', "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.strip('\n').split()
			if col[0] in genes.keys():
				period[col[0]] = col[1:]
		else:
			col = line.strip('\n').split('\t')
			periods = col[1:]

with open(source+'/Structure_ref.tdt', "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.strip('\n').split()
			if col[0] in genes.keys():
				structure[col[0]] = col[1:]
		else:
			col = line.strip('\n').split('\t')
			structures = col[1:]

with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split()
			gs = col[gene].split(':')
			gperiod = []
			gtissue = []
			gmaxp = []
			gmaxt = []
			for g in gs:
				gmaxperiod = "None"
				gmaxperiodvalue = 0
				gmaxtissue = "None"
				gmaxtissuevalue = 0
				if g == '-':
					gperiod.append('-')
					gtissue.append('-')
				else:
					p = []
					s = []
					for t in range(len(period[g])):
						if float(period[g][t]) > gmaxperiodvalue:
							gmaxperiodvalue = float(period[g][t])
							gmaxperiod = periods[t]
						if float(period[g][t]) >= float(thresh):
							p.append(t)
					for t in range(len(structure[g])):
						if float(structure[g][t]) > gmaxtissuevalue:
							gmaxtissuevalue = float(structure[g][t])
							gmaxtissue = structures[t]
						if float(structure[g][t]) >= float(thresh):
							s.append(t)
					ptis = []
					stis = []
					for m in p:
						ptis.append(periods[m])
					for m in s:
						stis.append(structures[m])
					if len(ptis) == 0:
						gperiod.append('-')
					else:
						gperiod.append(','.join(ptis))
					if len(stis) == 0:
						gtissue.append('-')
					else:
						gtissue.append(','.join(stis))
				gmaxp.append(gmaxperiod+"("+str(gmaxperiodvalue)+")")
				gmaxt.append(gmaxtissue+"("+str(gmaxtissuevalue)+")")
			ps = col[promoter].split(':')
			pperiod = []
			ptissue = []
			pmaxp = []
			pmaxt = []
			for g in ps:
				pmaxperiod = "None"
				pmaxperiodvalue = 0
				pmaxtissue = "None"
				pmaxtissuevalue = 0
				if g == '-':
					pperiod.append('-')
					ptissue.append('-')
				else:
					p = []
					s = []
					for t in range(len(period[g])):
						if float(period[g][t]) > pmaxperiodvalue:
							pmaxperiodvalue = float(period[g][t])
							pmaxperiod = periods[t]
						if float(period[g][t]) >= float(thresh):
							p.append(t)
					for t in range(len(structure[g])):
						if float(structure[g][t]) > pmaxtissuevalue:
							pmaxtissuevalue = float(structure[g][t])
							pmaxtissue = structures[t]
						if float(structure[g][t]) >= float(thresh):
							s.append(t)
					ptis = []
					stis = []
					for m in p:
						ptis.append(periods[m])
					for m in s:
						stis.append(structures[m])
					if len(ptis) == 0:
						pperiod.append('-')
					else:
						pperiod.append(','.join(ptis))
					if len(stis) == 0:
						ptissue.append('-')
					else:
						ptissue.append(','.join(stis))
				pmaxp.append(pmaxperiod+"("+str(pmaxperiodvalue)+")")
				pmaxt.append(pmaxtissue+"("+str(pmaxtissuevalue)+")")
			print line.strip('\n')+"\t"+':'.join(gmaxp)+"\t"+':'.join(gperiod).strip('\n')+"\t"+':'.join(gmaxt)+"\t"+':'.join(gtissue).strip('\n')+"\t"+':'.join(pmaxp)+"\t"+':'.join(pperiod).strip('\n')+"\t"+':'.join(pmaxt)+"\t"+':'.join(ptissue).strip('\n')
		else:
			print line.strip('\n')+"\tGene_max_brainspan_period\tGene_brainspan_periods(at_least_"+str(thresh)+")\tGene_max_brainspan_tissue\tGene_brainspan_tissues(at_least_"+str(thresh)+")\tPromoter_max_brainspan_period\tPromoter_brainspan_periods(at_least_"+str(thresh)+")\tPromoter_max_brainspan_tissue(at_least_"+str(thresh)+")\tPromoter_brainspan_tissues"


'''
#!/bin/bash
ann_folder=$2
thresh=$3
function join_by { local IFS="$1"; shift; echo "$*"; }

hed=$(head -n 1 $1)
while read line ; do
	if [[ $line =~ ^# ]]; then
		g_field=$(echo $line | sed 's/ /\n/g' | grep -nx gene | cut -d: -f1)
		p_field=$(echo $line | sed 's/ /\n/g' | grep -nx promoter | cut -d: -f1)
		echo -e "${line}\tGene_max_brainspan_period\tGene_brainspan_periods\tGene_max_brainspan_tissue\tGene_brainspan_tissues\tPromoter_max_brainspan_period\tPromoter_brainspan_periods\tPromoter_max_brainspan_tissue\tPromoter_brainspan_tissues"
	else
		g_per=''
		g_M_per=''
		g_tis=''
		g_M_tis=''
		p_per=''
		p_M_per=''
		p_tis=''
		p_M_tis=''
		genes=$(echo $line | cut -d' ' -f$g_field | sed 's/:/\n/g')
		promoters=$(echo $line | cut -d' ' -f$p_field | sed 's/:/\n/g')
		for gene in $genes; do
			if [ $gene != '-' ] && [ $(grep ${gene} ${ann_folder}/Period_ref.tdt | head -n 1 | wc -l | cut -d: -f1) = '1' ]; then
				per=''
				tis=''
				per_line=$(awk -v Gen=${gene} '{if ($1==Gen) print $0}' ${ann_folder}/Period_ref.tdt)
				tis_line=$(awk -v Gen=${gene} '{if ($1==Gen) print $0}' ${ann_folder}/Structure_ref.tdt)
				max_per=$(echo -e $per_line | sed 's/ /\n/g' | grep -nx $(echo -e $per_line | sed 's/ /\n/g'| sort -r -n | head -n 1) | cut -d: -f1)
				max_tis=$(echo -e $tis_line | sed 's/ /\n/g' | grep -nx $(echo -e $tis_line | sed 's/ /\n/g'| sort -r -n | head -n 1) | cut -d: -f1)
				if (( $(echo "$(echo -e $per_line | cut -d' ' -f$max_per) > $thresh" | bc -l ) )); then
					M_per=$(echo "$(head -n 1 ${ann_folder}/Period_ref.tdt | cut -f$max_per)($(echo -e $per_line | cut -d' ' -f$max_per))")
					for col in $(echo $per_line | cut -d' ' -f2-); do
						if (( $(echo "$col > $thresh" | bc -l) )); then
							per=$(echo "${per},$(head -n 1 ${ann_folder}/Period_ref.tdt | cut -f$(echo -e $per_line | sed 's/ /\n/g' | grep -nx ${col} | cut -d: -f1))(${col})")
						fi
					done
				else
					M_per='-'
					per=',-'
				fi
				per=${per#?}
				if (( $(echo "$(echo -e $tis_line | cut -d' ' -f$max_tis) > $thresh" | bc -l ) )); then
					M_tis=$(echo "$(head -n 1 ${ann_folder}/Structure_ref.tdt | cut -f$max_tis)($(echo -e $tis_line | cut -d' ' -f$max_tis))")
					for col in $(echo $tis_line | cut -d' ' -f2-); do
						if (( $(echo "$col > $thresh" | bc -l) )); then
							tis=$(echo "${tis},$(head -n 1 ${ann_folder}/Structure_ref.tdt | cut -f$(echo -e $tis_line | sed 's/ /\n/g' | grep -nx ${col} | cut -d: -f1))(${col})")
						fi
					done
				else
					M_tis='-'
					tis=',-'
				fi
				tis=${tis#?}
				g_per=$(echo "${g_per}:${per}")
				g_M_per=$(echo "${g_M_per}:${M_per}")
				g_tis=$(echo "${g_tis}:${tis}")
				g_M_tis=$(echo "${g_M_tis}:${M_tis}")
			else
				g_per=$(echo "${g_per}:-")
				g_M_per=$(echo "${g_M_per}:-")
				g_tis=$(echo "${g_tis}:-")
				g_M_tis=$(echo "${g_M_tis}:-")
			fi
		done
		for gene in $promoters; do
			if [ $gene != '-' ] && [ $(grep ${gene} ${ann_folder}/Period_ref.tdt | head -n 1 | wc -l | cut -d: -f1) = '1' ]; then
				per=''
				tis=''
				per_line=$(awk -v Gen=${gene} '{if ($1==Gen) print $0}' ${ann_folder}/Period_ref.tdt)
				tis_line=$(awk -v Gen=${gene} '{if ($1==Gen) print $0}' ${ann_folder}/Structure_ref.tdt)
				max_per=$(echo -e $per_line | sed 's/ /\n/g' | grep -nx $(echo -e $per_line | sed 's/ /\n/g'| sort -r -n | head -n 1) | cut -d: -f1)
				max_tis=$(echo -e $tis_line | sed 's/ /\n/g' | grep -nx $(echo -e $tis_line | sed 's/ /\n/g'| sort -r -n | head -n 1) | cut -d: -f1)
				if (( $(echo "$(echo -e $per_line | cut -d' ' -f$max_per) > $thresh" | bc -l ) )); then
					M_per=$(echo "$(head -n 1 ${ann_folder}/Period_ref.tdt | cut -f$max_per)($(echo -e $per_line | cut -d' ' -f$max_per))")
					for col in $(echo $per_line | cut -d' ' -f2-); do
						if (( $(echo "$col > $thresh" | bc -l) )); then
							per=$(echo "${per},$(head -n 1 ${ann_folder}/Period_ref.tdt | cut -f$(echo -e $per_line | sed 's/ /\n/g' | grep -nx ${col} | cut -d: -f1))(${col})")
						fi
					done
				else
					M_per='-'
					per=',-'
				fi
				per=${per#?}
				if (( $(echo "$(echo -e $tis_line | cut -d' ' -f$max_tis) > $thresh" | bc -l ) )); then
					M_tis=$(echo "$(head -n 1 ${ann_folder}/Structure_ref.tdt | cut -f$max_tis)($(echo -e $tis_line | cut -d' ' -f$max_tis))")
					for col in $(echo $tis_line | cut -d' ' -f2-); do
						if (( $(echo "$col > $thresh" | bc -l) )); then
							tis=$(echo "${tis},$(head -n 1 ${ann_folder}/Structure_ref.tdt | cut -f$(echo -e $tis_line | sed 's/ /\n/g' | grep -nx ${col} | cut -d: -f1))(${col})")
						fi
					done
				else
					M_tis='-'
					tis=',-'
				fi
				tis=${tis#?}
				p_per=$(echo "${p_per}:${per}")
				p_M_per=$(echo "${p_M_per}:${M_per}")
				p_tis=$(echo "${p_tis}:${tis}")
				p_M_tis=$(echo "${p_M_tis}:${M_tis}")
			else
				p_per=$(echo "${p_per}:-")
				p_M_per=$(echo "${p_M_per}:-")
				p_tis=$(echo "${p_tis}:-")
				p_M_tis=$(echo "${p_M_tis}:-")
			fi
		done
		g_per=${g_per#?}
		g_M_per=${g_M_per#?}
		g_tis=${g_tis#?}
		g_M_tis=${g_M_tis#?}
		p_per=${p_per#?}
		p_M_per=${p_M_per#?}
		p_tis=${p_tis#?}
		p_M_tis=${p_M_tis#?}
		echo -e "${line}\t$g_M_per\t$g_per\t$g_M_tis\t$g_tis\t$p_M_per\t$p_per\t$p_M_tis\t$p_tis"
	fi
done < $1
'''