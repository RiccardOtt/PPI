import sys
import numpy as np





def mut_converter(mutation):
	new_list = []
	final_list = []
	real_final_list = []
	convertion = ''

	d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	for k,v in d.items():
		if mutation[0] == v:
			mutation[0] = k
			new_list.append(mutation)

	for i in new_list:
		for k,v in d.items():
			i = list(i)
			if i[-1] == v:
				i[-1] = k
				final_list.append(i)

	for a in final_list:
#		print(a)
		real_final_list.append(mutation[0]+' '+a[1]+' '+a[2]+''+''.join(a[3:-1])+' '+a[-1])

	for f in real_final_list:
		convertion += f

	return convertion





def correlation(external_pred,mypred):

	my_ddg = []
	for h in mypred:
		h = h.split('   ')
		my_ddg.append(h[0])

	check = []
	for j in mypred:
		j = j.split('   ')
		print(j)
		check.append(j[0])


	common_ddg = []
	for i in external_pred:
		i = i.split('   ')
		if i[0] in check:
#			print(i)
			common_ddg.append(i[1])


#	extarr = np.array(common_ddg)
#	intarr = np.array(my_ddg)
#	corr = np.corrcoef(extarr,intarr)
#	print(corr)









if __name__ == '__main__':
	pred = sys.argv[1]
	my_pred = sys.argv[2]

	pred_conv = []
	my_predct = []


	with open(pred) as ext_pred:
		for lines in ext_pred:
			lines = lines.rstrip().split(' ')
			mut = lines[1]
			mut = list(mut)
			mut_conv = mut_converter(mut)
			pred_conv.append(lines[0]+' '+mut_conv+'   '+'DDG_pred = '+lines[2])


	with open(my_pred) as my_predictions:
		for lines in my_predictions:
			lines = lines.rstrip()
			my_predct.append(lines)

	correlation(pred_conv,my_predct)
