ann_file=$2
function join_by { local IFS="$1"; shift; echo "$*"; }

while read line ; do
	if [[ $line =~ ^# ]]; then
		g_field=$(echo $line | sed 's/ /\n/g' | grep -nx gene | cut -d: -f1)
		p_field=$(echo $line | sed 's/ /\n/g' | grep -nx promoter | cut -d: -f1)
		echo -e "${line}\tGene_HPO_phenotypes\tPromoter_HPO_phenotypes"
	else
		g_ann=''
		p_ann=''
		genes=$(echo $line | cut -d' ' -f$g_field | sed 's/:/\n/g')
		promoters=$(echo $line | cut -d' ' -f$p_field | sed 's/:/\n/g')
		for gene in $genes; do
			if [ $gene != '-' ]; then
				ann=$( join_by / $(fgrep -v "#" $ann_file | awk -v gen=$gene '{if ($2 == gen) print $0}' | sed 's/[ -]/_/g' | cut -f3))
				if [ -z $ann ]; then
					ann='-'
				fi
				g_ann=$(echo -e "${g_ann}:${ann}")
			fi
		done
		for promoter in $promoters; do
			if [ $promoter != '-' ]; then
				ann=$( join_by / $(fgrep -v "#" $ann_file | awk -v gen=$promoter '{if ($2 == gen) print $0}' | sed 's/[ -]/_/g' | cut -f3))
				if [ -z $ann ]; then
					ann='-'
				fi
				p_ann=$(echo -e "${p_ann}:${ann}")
			fi
		done
		g_ann=${g_ann#?}
		p_ann=${p_ann#?}
		if [ -z $g_ann ]; then
			g_ann='-'
		fi
		if [ -z $p_ann ]; then
			p_ann='-'
		fi
		echo -e "${line}\t${g_ann}\t${p_ann}"
	fi
done < $1
