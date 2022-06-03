from sys import argv
import sys
import os
script, a, b, c = argv

def Average(lst):
	return sum(lst) / len(lst)

genes=[]
with open(a, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split('"')
			genes.append(col[3])

samps={}
with open(b, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split(',')
			try:
				samps[col[2]]
			except KeyError:
				samps[col[2]] = []
			samps[col[2]].append(col[0])
head=["Gene"]
for samp in samps:
	head.append(str(samp))
	#head.append(str(samp+"_min"))
	#head.append(str(samp+"_max"))
print "\t".join(head)
cnt = 0
with open(c, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			ret = [str(genes[cnt])]
			col = line.split(',')
			for samp in samps:
				nums=[]
				for x in samps[samp]:
					nums.append(float(col[int(x)]))
				nums.sort()
				#ret.append(str(Average(nums)))
				ret.append(str(nums[0]))
				#ret.append(str(nums[len(nums)-1]))
			print "\t".join(ret)
			cnt+=1


'''
#!/bin/bash

samples=$(cut -d',' -f3 $1 | sort -u)
hed='#line'
for s in $samples;do
	hed=$(echo -e "$hed\t${s}_min\t${s}_av\t${s}_max")
done
echo $hed
while read line;do
	ret=$( head -n $(( $(echo $line | cut -d',' -f1)+1 )) $3 | tail -n 1 | cut -d'"' -f4)
	for s in $samples;do
		vals=''
		for l in $(awk -v q="'" --field-separator ',' -v sam=$s '{if ($3==sam) print $1+1}' $1); do
			vals=$(echo "$vals|$(echo $line | cut -d',' -f$l)")
		done
		vals=${vals#?}
		max=$(echo $vals | sed 's/|/\n/g' | sort -n -r | head -n 1)
		min=$(echo $vals | sed 's/|/\n/g' | sort -n | head -n 1)
		av=$(echo $vals | sed 's/|/\n/g' | awk '{ total += $1; count++ } END { print total/count }')
		ret=$(echo -e "$ret\t$max\t$min\t$av")
	done
	echo $ret
done <$2
'''
