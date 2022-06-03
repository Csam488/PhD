from sys import argv
import sys
import os
script, a = argv

def Average(lst):
	return sum(lst) / len(lst)

genes=[]
ref = {}
ref['A1C']= []
with open(a+'A1C.tdt', "r") as fl:
	for cnt, line in enumerate(fl):
		if (cnt > 0):
			col = line.split()
			genes.append(col.pop(0))
			num = []
			for c in col:
				num.append(float(c))
			ref['A1C'].append(Average(num))
files=['AMY','CBC','CB','CGE','DFC','DTH','HIP','IPC','ITC','LGE','M1C-S1C','M1C','MD','MFC','MGE','Ocx','OFC','PCx','S1C','STC','STR','TCx','URL','V1C','VFC']
for f in files:
	ref[f]= []
	with open(a+f+'.tdt', "r") as fl:
		for cnt, line in enumerate(fl):
			if (cnt > 0):
				col = line.split()
				col.pop(0)
				num = []
				for c in col:
					num.append(float(c))
				ref[f].append(Average(num))

head = '#Gene\tA1C'
for f in files:
	head=head+'\t'+f
print head

for x in range(0,len(genes)):
	ret = genes[x]+'\t'+str(ref['A1C'][x])
	for f in files:
		ret=ret+'\t'+str(ref[f][x])
	print ret
	#sys.exit()
