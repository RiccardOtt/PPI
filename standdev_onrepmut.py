import sys
import os
import numpy as np
from collections import defaultdict


def stat(file):
	aff_wt = []
	aff_mut = []
	rep_mut = []
	rep_id = []
	mut = []
	d_wt = {}
	d_mut = {}
	std_wt = []
	std_mut = []

	for lines in file:
		lines = lines.rstrip().split(';')
		if ',' not in lines[1]:
			rep_mut.append(lines[1])
			rep_id.append(lines[0])
			aff_wt.append(lines[8])
			aff_mut.append(lines[6])


	for i in range(len(rep_mut)):
		d_wt[rep_mut[i]] = []


	for keys,values in d_wt.items():
		for i in range(len(rep_mut)):
			if keys == rep_mut[i]:
				d_wt[rep_mut[i]].append(aff_wt[i])


	for k,v in d_wt.items():
#		print(k)
		arr = np.array(v).astype(np.float)
		std = np.std(arr)
		std_wt.append(std)


	for i in range(len(rep_mut)):
		d_mut[rep_mut[i]] = []


	for keys,values in d_mut.items():
		for i in range(len(rep_mut)):
			if keys == rep_mut[i]:
				d_mut[rep_mut[i]].append(aff_mut[i])


	for k,v in d_mut.items():
		arr = np.array(v).astype(np.float)
		std = np.std(arr)
		std_mut.append(std)

	print(std_wt)
#	print(std_mut)

	array_wt = np.std(std_wt)
#	print(array_wt)


	array_mut = np.std(std_mut)
#	print(array_mut)



if __name__ == '__main__':
	file = sys.argv[1]
	parse = open(file)
	stat(parse)
