import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import matplotlib.ticker as mtick
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set_theme()
import pandas as pd


def heatmap(statmatrix):
	x_axis_labels = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	y_axis_labels = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	arr = np.loadtxt(statmatrix)
#	ax = sns.heatmap(arr, linewidth=0.5)
#	plt.show()
	sns.heatmap(arr,cmap='coolwarm', xticklabels=x_axis_labels, yticklabels=y_axis_labels)
	plt.show()




def histogram(interactres,surfres):
	res = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	surfr = 0
	interesid = 0
	listsurf = []
	listint = []


	dsurf = { i : 0 for i in res}
	for lines in surfres:
		surfr += 1
		lines = lines.rstrip()
		if lines in res:
			dsurf[lines] += 1
	for k,v in dsurf.items():
		dsurf[k] = (v/surfr)*100


	dinter = { i : 0 for i in res}
	for lines in interactres:
		interesid += 1
		lines = lines.rstrip().split(' ')
		r1 = lines[1]
		r2 = lines[7]
		if r1 in res:
			if r2 in res:
				dinter[r1] += 1
				dinter[r2] += 1
	for k,v in dinter.items():
		dinter[k] = (v/interesid)*100


	for k,v in dsurf.items():
		listsurf.append(v)

	for k,v in dinter.items():
		listint.append(v)


	width = 0.17
	x = np.arange(len(res))
	x1 = [a + width for a in x]
	x2 = [a + width for a in x1]
#	x3 = [a + width for a in x2]

#	plt.bar(x, H_norm, width, color='yellow', label='H')
	plt.bar(x1, listsurf, width, color='blue', label='surface residues')
	plt.bar(x2, listint, width, color='red', label='interacting residues')
#	plt.bar(x3, overall_norm, width, color='lightsteelblue', label='Overall')

	plt.ylabel('Composition(%)')
	plt.xlabel('Aminoacid')
	plt.yticks()
#	plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
	plt.xticks(x,res,rotation='vertical')
	plt.legend(loc="upper right")
#	plt.show()





if __name__ == '__main__':
	matrix = sys.argv[1]
	intres = sys.argv[2]
	surfaceres = sys.argv[3]
	with open(matrix) as m, open(intres) as interacres, open(surfaceres) as surfres:
		heatmap(m)
#		histogram(interacres,surfres)
