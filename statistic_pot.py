import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def interacting_dict(interacting_residue):

	interacting_set = set()
	interacting_dict = {}
	interacting_res = []
	residues_dict = {}
	tot_int_pairs = 0

	for lines in interacting_residue:
		if 'UNK' not in lines:
			if 'SEP' not in lines:
				if 'TPO' not in lines:
					if 'GVE' not in lines:
						if 'F2F' not in lines:
							tot_int_pairs += 1
							lines = lines.rstrip().split(' ')
							interacting_set.add(lines[1]+' '+lines[10])
							interacting_res.append(lines[1]+' '+lines[10])

	interacting_dict = {k:0 for k in interacting_set}

	for int in interacting_res:
		int = int.split(' ')
		for res in int:
			if res != 'UNK':
				if res != 'SEP':
					if res != 'TPO':
						if res != 'GVE' and res != 'F2F':
							if res not in residues_dict:
								residues_dict[res] = 0


#	print(interacting_res)
	return interacting_dict,interacting_res,tot_int_pairs,residues_dict




def Matrix(contacts):

	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	M = []

	for pair in contacts:
		pair = pair.split(' ')
		for res in pair:
			if res not in amino:
				amino.append(res)

	for i in range(len(amino)+1):
		M.append([])
		for j in range(len(amino)+1):
			M[i].append(0)

#	M = np.array(M)
#	M = pd.DataFrame(M, columns=list(amino), index=list(amino))

	for i in range(len(amino)):
		M[i+1][0] = amino[i]
		M[0][i+1] = amino[i]


#	for k in M:
#		print(k)

	return M





def statistic_pot(contacts,cont_dict,res_dict,allresidues,tot_inter_pairs,matrix):
	mol_fract = {}

	Tot_res = 0

	mol_fract = {k:0 for k in res_dict}

	Stat_Pot = set()

	for i in contacts:
		for keys in cont_dict:
			if i == keys:
				cont_dict[i] += 1

	for lines in allresidues:
		residue = lines.rstrip()
		Tot_res += 1
		for k,v in res_dict.items():
			if residue == k:
				res_dict[k] +=1

#	print(res_dict)

	for k,v in res_dict.items():
#		print(res_dict)
		for K,V in mol_fract.items():
			mol_fract[k] = v/Tot_res


	for contact in contacts:
		contact = contact.split(' ')
#		print(cont_dict[contact[0]+' '+contact[1]])
		P = math.log((cont_dict[contact[0]+' '+contact[1]])/(mol_fract[contact[0]]*mol_fract[contact[1]]*tot_inter_pairs))
#		print(contact[0],contact[1],P)
		Stat_Pot.add(contact[0]+' '+contact[1]+' '+str(P))



	for i in range(len(matrix)):
		for j in range(len(matrix)):
			for r in Stat_Pot:
				r = r.split(' ')
				if matrix[i][0] == r[0]:
					if matrix[0][j] == r[1]:
						matrix[i][j] = round(float(r[2]),1)

#	for i in Stat_Pot:
#		print(i)
#	for h in matrix:
#		print(h)

	pretty_matrix = pd.DataFrame(matrix)
#	pretty_matrix.to_csv('my_matrix_7A_relacc.csv',sep='\t')

	matrix = np.array(matrix)
	pretty_matrix.to_csv('matrix_intres_6A_allres.txt', sep='\t')

	print(pretty_matrix)





if __name__ == '__main__':
	interacting_res = sys.argv[1]
	all_res = sys.argv[2]
	int_res = open(interacting_res,'r')
	allresidues = open(all_res,'r')

	inter_dict,inter_res,tot_interact_pairs,residues_dict = interacting_dict(int_res)
	empty_M = Matrix(inter_res)
	statistic_pot(inter_res,inter_dict,residues_dict,allresidues,tot_interact_pairs,empty_M)
