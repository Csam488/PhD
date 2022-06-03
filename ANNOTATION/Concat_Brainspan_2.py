from sys import argv
import sys
import os
script, a = argv

def Average(lst):
	return sum(lst) / len(lst)

genes=[]
ref = {}
ref['Br_2A']= []
cols=[]
with open(a+'Br_2A.tdt', "r") as fl:
	for cnt, line in enumerate(fl):
		if (cnt > 0):
			col = line.split()
			genes.append(col[0])
			num = []
			for c in cols:
				num.append(float(col[c]))
			ref['Br_2A'].append(Average(num))
		else:
			col = line.split()
			for c in range(0,len(col)):
				if '_average' in col[c]:
					cols.append(c)
files=['Br_2B','Br_3A','Br_3B','Br_4','Br_5','Br_6','Br_7','Br_8','Br_9','Br_10','Br_11']
for f in files:
	ref[f]= []
	cols=[]
	with open(a+f+'.tdt', "r") as fl:
		for cnt, line in enumerate(fl):
			if (cnt > 0):
				col = line.split()
				num = []
				for c in cols:
					num.append(float(col[c]))
				ref[f].append(Average(num))
			else:
				col = line.split()
				for c in range(0,len(col)):
					if '_average' in col[c]:
						cols.append(c)

head = '#Gene\tBr_2A'
for f in files:
	head=head+'\t'+f
print head

for x in range(0,len(genes)):
	ret = genes[x]+'\t'+str(ref['Br_2A'][x])
	for f in files:
		ret=ret+'\t'+str(ref[f][x])
	print ret
	#sys.exit()
