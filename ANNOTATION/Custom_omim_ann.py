import sys
from sys import argv
import re
script, input, source = argv

genes = {}

with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.strip('\n').split()
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

for g in genes.keys():
	genes[g] = []

with open(source, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.strip('\n').split('\t')
			if col[1] in genes.keys():
				genes[col[1]].append(col[2].replace(',',' ').replace(':',' '))

for g in genes.keys():
	if len(genes[g]) == 0:
		genes[g] = '-'
	else:
		genes[g] = ','.join(genes[g])
genes['-'] = '-'

with open(input, "r") as fl:
	for cnt, line in enumerate(fl):
		if line[0] != '#':
			col = line.split()
			gs = col[gene].split(':')
			gout = []
			for g in gs:
				gout.append(genes[g])
			ps = col[promoter].split(':')
			pout = []
			for p in ps:
				pout.append(genes[p])
			print line.strip('\n')+'\t'+':'.join(gout)+'\t'+':'.join(pout)
		else:
			print line.strip('\n')+"\tGene_HPO_phenotypes\tPromoter_HPO_phenotypes"
