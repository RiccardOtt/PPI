import sys
import numpy as np
import pandas as pd
import Bio.PDB
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)



def mut_converter(skempi_line):
	new_list = []
	final_list = []
	real_final_list = []
	convertion = ''

	d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	for k,v in d.items():
		mutation = list(skempi_line[2])
		if mutation[1] == v:
			mutation[1] = k
			new_list.append(mutation)

	for i in new_list:
		for k,v in d.items():
			i = list(i)
			if i[-1] == v:
				i[-1] = k
				final_list.append(i)

	for a in final_list:
		real_final_list.append(skempi_line[0]+' '+a[1]+' '+''.join(a[3:-1])+' '+a[-1])

	for f in real_final_list:
		convertion += f

	return convertion




def interface(id):

	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	contacts = []
	parser = PDBParser()

	structure = parser.get_structure('PHA-L', '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/'+id+'.pdb')


	model = structure[0]
	for chain1 in model:
		for chain2 in model:
			if chain1 != chain2:
				for residue1 in chain1:
					for residue2 in chain2:
						tags1 = residue1.id
						tags2 = residue2.id
						if tags1[0] != " " or tags2[0] != " ":
							pass
						else:
							try:
								distance = residue1['CA'] - residue2['CA']
							except KeyError: ## no CA atom, e.g. for H_NAG
								continue
							if distance < 8:
								cont = str(residue1)+' '+str(residue2)
								if cont not in contacts:
									contacts.append(str(residue1)+' '+str(chain1)+' '+str(residue2)+' '+str(chain2))

	for i in contacts:
		print(i)


	return contacts





def filt_relsurf_acc(interf_contacts,id):

	cont_filt = []
	for lines in interf_contacts:
		lines = lines.split(' ')

		try:
			chainA = list(lines[8])[3]
			resA = lines[4].split('=')[1]
		except:	continue

		try:
			chainB = list(lines[17])[3]
			resB = lines[13].split('=')[1]
		except:	continue

		try:
			dssp_file_complex = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/complex/'+id+'.dssp','r')
			dssp_file_singlchainA = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/only_chain/'+id+'_'+chainA+'.dssp','r')
			dssp_file_singlchainB = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/only_chain/'+id+'_'+chainB+'.dssp','r')
		except:	continue
		else:
			rel_acc = []
			for i, line in enumerate(dssp_file_complex):
				if i > 27:
					line = line.rstrip()
					n_line = line.split(' ')
					n_line = list(filter(None, n_line))

					if resA == n_line[1]:
						if chainA == n_line[2]:
							rel_acc.append(line[35:38])

					if resB == n_line[1]:
						if chainB == n_line[2]:
							rel_acc.append(line[35:38])


			for i, line in enumerate(dssp_file_singlchainA):
				if i > 27:
					line = line.rstrip()
					n_line = line.split(' ')
					n_line = list(filter(None, n_line))

					if resA == n_line[1]:
						if chainA == n_line[2]:
							rel_acc.append(line[35:38])


			for i, line in enumerate(dssp_file_singlchainB):
				if i > 27:
					line = line.rstrip()
					n_line = line.split(' ')
					n_line = list(filter(None, n_line))

					if resB == n_line[1]:
						if chainB == n_line[2]:
							rel_acc.append(line[35:38])
		try:

			if abs(float(rel_acc[0]))-float(float(rel_acc[1])) < abs(float(rel_acc[2])-float(rel_acc[3])):
				cont_filt.append(lines[0:16])
		except:
			pass


#	for i in cont_filt:
#		print(i)
	return cont_filt






def intracontacts(id):

	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	intracontacts = []
	parser = PDBParser()

	structure = parser.get_structure('PHA-L', '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/'+id+'.pdb')

	model = structure[0]
	for chain1 in model:
		for chain2 in model:
			for residue1 in chain1:
				for residue2 in chain2:
				# compute distance between CA atoms
					try:
						distance = residue1['CA'] - residue2['CA']
					except KeyError:
       				## no CA atom, e.g. for H_NAG
						continue
					if distance < 6:
						cont = str(residue1)+' '+str(residue2)
						if cont not in intracontacts:
							intracontacts.append(cont)


#	for i in intracontacts:
#		print(i)
	return intracontacts





def prediction(matrix1,matrix2,intercont,intracont):
	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	scoringdict1 = {}
	scoringdict2 = {}
	intradict = {}
	interdict = {}




	for lines in intracont:
		lines = lines.split(' ')
		intrares = lines[1]+' '+lines[8]
		if intrares not in intradict:
			intradict[intrares] = 1
		else:
			intradict[intrares] += 1

	for lines in intercont:
		interres = (lines[1]+' '+lines[10])
#		print(interres)
		if interres not in interdict:
			interdict[interres] = 1
		else:
			interdict[interres] += 1


	df1 = pd.DataFrame(matrix1,amino,amino)
	df2 = pd.DataFrame(matrix2,amino,amino)
	l = []
	for i in amino:
		for g in amino:
			p = i+' '+g
			scoringdict1[p] = df1[i][g]
			scoringdict2[p] = df2[i][g]

#	print(scoringdict2)

	interscore = 0
	intrascore = 0
	for k,v in interdict.items():
		for K,V in scoringdict1.items():
			if k == K:
				interscore += (v*V)/len(interdict)

	for k1,v1 in intradict.items():
		for K1,V1 in scoringdict2.items():
			if k1 == K1:
				intrascore += (v1*V1)/len(intradict)


#	print('interscore = ',interscore)
#	print('intrascore = ',intrascore)





if __name__ == '__main__':
	matrix1 = sys.argv[1]
	matrix2 = sys.argv[2]
	skempi_single = sys.argv[3]
	skem = open(skempi_single)
	id_list = []
	mut_list = []


	matrix1 = np.loadtxt(matrix1)
	matrix2 = np.loadtxt(matrix2)
#	matrix2 = open(matrix2)
#	matrix2 = matrix2.read()
#	matrix2 = np.fromstring(matrix2, sep='\n')
#	matrix2 = list(matrix2)
#	print(matrix1)


	for lines in skem:
		lines = lines.split(',')
		mutation = mut_converter(lines)
		mut_list.append(mutation)
		pdb_id = lines[0]
		ids = pdb_id.split('_')
		id = ids[0]
		if id not in id_list:
			id_list.append(id)


	for pdb_id in id_list:
		intercontacts = interface(pdb_id)
		filt_relacc_intercontacts = filt_relsurf_acc(intercontacts,pdb_id)
		intracont = intracontacts(pdb_id)
		prediction(matrix1,matrix2,filt_relacc_intercontacts,intracont)
		break
