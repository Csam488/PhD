ann_folder=$2
function join_by { local IFS="$1"; shift; echo "$*"; }
#module load BEDTools/2.28.0-gimkl-2018b

while read chr st sp oth ; do
	if [[ $chr =~ ^# ]]; then
		echo -e "${chr}\t${st}\t${sp}\t${oth}\tReg_element\tReg_element_loc"
	else
		H=$( join_by , $( echo -e "${chr}\t${st}\t${sp}\n" | bedtools intersect -b - -a ${ann_folder}/Reg_ele.bed | awk -v a=':' -v b='-' '{print $4}'))
		M=$( join_by , $( echo -e "${chr}\t${st}\t${sp}\n" | bedtools intersect -b - -a ${ann_folder}/Reg_ele.bed | awk -v a=':' -v b='-' '{print $1a$2b$3}'))
		if test -z "$M" ;then
			M='-'
		elif [ $(IFS=','; set -f; set -- $M; echo $#) -gt 10 ];then
			M='More_than_10'
		fi
		if test -z "$H" ;then
			H='-'
		elif [ $(IFS=','; set -f; set -- $H; echo $#) -gt 10 ];then
			H='More_than_10'
		fi
		echo -e "${chr}\t${st}\t${sp}\t${oth}\t${H}\t${M}"
	fi
done < $1
