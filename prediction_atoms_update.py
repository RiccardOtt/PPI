import sys
import os
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
#		print(a)
		real_final_list.append(skempi_line[0]+' '+a[1]+' '+a[2]+' '+''.join(a[3:-1])+' '+a[-1])

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
			for residue1 in chain1:
				tags1 = residue1.id
#				print(residue1)
				for residue2 in chain2:
					tags2 = residue2.id
					if chain1 != chain2:
						if tags1[0] != " " or tags2[0] != " ":
							pass
						else:
#						try:
							atoms1 = residue1.get_atoms()
							atoms2 = residue2.get_atoms()
							for atom1 in atoms1:
								for atom2 in atoms2:
									distance = atom1 -atom2
#									print(distance)
#							distance = residue1['CA'] - residue2['CA']
#									except KeyError:
         	       				## no CA atom, e.g. for H_NAG
#										continue
									if distance < int(8):
#										print(residue1,residue2)
										contact = (str(residue1)+' '+str(chain1)+' '+str(residue2)+' '+str(chain2))
										if contact not in contacts:
											contacts.append(str(residue1)+' '+str(chain1)+' '+str(residue2)+' '+str(chain2))

#	for i in contacts:
#		print(i)


	return contacts





def filt_relsurf_acc(interf_contacts,id):

	cont_filt = []
	for lines in interf_contacts:
		lines = lines.split(' ')
#		print(lines)

		try:
			chainA = list(lines[8])[3]
			resA = lines[4].split('=')[1]
#			print(resA,chainA)
		except:	continue

		try:
			chainB = list(lines[17])[3]
			resB = lines[13].split('=')[1]
#			print(resB,chainB)
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
			if chain1 == chain2:
				for residue1 in chain1:
					tags1 = residue1.id
					for residue2 in chain2:
						tags2 = residue2.id
						if tags1[0] != " " or tags2[0] != " ":
				  	  ###	The residue is a heteroatom
							pass
						else:
							try:
				# compute distance between CA atoms
								distance = residue1['CA'] - residue2['CA']
      				## no CA atom, e.g. for H_NAG
							except KeyError:
								continue

							if distance < 8:
								cont = str(residue1)+' '+str(chain1)+' '+str(residue2)+' '+str(chain2)
								if cont not in intracontacts:
									intracontacts.append(cont)


#	for i in intracontacts:
#		print(i)
	return intracontacts






def prediction(matrix1,matrix2,intercont,intracont,mut):
	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	amino2 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

	wt_intracont = []
	wt_intercont = []
	mut_intracont = []
	mut_intercont = []

	scoringdict1 = {}
	scoringdict2 = {}

	wt_intradict = {}
	wt_interdict = {}

	mut_scoringdict1 = {}
	mut_scoringdict2 = {}
	mut_intradict = {}
	mut_interdict = {}



	for lines in intracont:
		lines = lines.split(' ')
		intrares = lines[1]+' '+lines[4].split('=')[1]+' '+lines[8].split('=')[1].replace('>','')+' '+lines[10]
		wt_intracont.append(intrares.split(' ')[0]+' '+intrares.split(' ')[3])
		if mut.split(' ')[3] == intrares.split(' ')[1] and mut.split(' ')[1] == intrares.split(' ')[0]:
			intrares = intrares.replace(intrares.split(' ')[0],mut.split(' ')[4])
		mut_intracont.append(intrares.split(' ')[0]+' '+intrares.split(' ')[3])

	for lines in intercont:
		interres = (lines[1]+' '+lines[10])
		interres = lines[1]+' '+lines[4].split('=')[1]+' '+lines[8].split('=')[1].replace('>','')+' '+lines[10]
		wt_intercont.append(interres.split(' ')[0]+' '+interres.split(' ')[3])
		if mut.split(' ')[3] == interres.split(' ')[1] and mut.split(' ')[1] == interres.split(' ')[0]:
			interres = interres.replace(interres.split(' ')[0],mut.split(' ')[4])
		mut_intercont.append(interres.split(' ')[0]+' '+interres.split(' ')[3])


	for wt in wt_intracont:
		if wt not in wt_intradict:
			wt_intradict[wt] = 1
		else:
			wt_intradict[wt] += 1
	for mt in mut_intracont:
		if mt not in mut_intradict:
			mut_intradict[mt] = 1
		else:
			mut_intradict[mt] += 1


	for wt in wt_intercont:
		if wt not in wt_interdict:
			wt_interdict[wt] = 1
		else:
			wt_interdict[wt] += 1
	for mt in mut_intercont:
		if mt not in mut_interdict:
			mut_interdict[mt] = 1
		else:
			mut_interdict[mt] += 1



	df1 = pd.DataFrame(matrix1,amino2,amino2)
	df2 = pd.DataFrame(matrix2,amino,amino)
	l = []

	for i in amino:
		for g in amino:
			p = i+' '+g
#			scoringdict1[p] = df1[i][g]
			scoringdict2[p] = df2[i][g]

	for h in amino2:
		for s in amino2:
			P = h+' '+s
			scoringdict1[P] = df1[h][s]
#	print(scoringdict1)

	tot1 = 0
	tot2 = 0
	tot3 = 0
	tot4 = 0
	wt_intrascore = 0
	wt_interscore = 0
	mut_intrascore = 0
	mut_interscore = 0

	for k,v in wt_intradict.items():
		tot1 += v
		for K,V in scoringdict1.items():
			if k == K:
				wt_intrascore += (v*V) #/len(wt_intradict)

	for k1,v1 in wt_interdict.items():
		tot2 += v
		for K1,V1 in scoringdict2.items():
			if k1 == K1:
				wt_interscore += (v1*V1) #/len(wt_interdict)

	for k,v in mut_intradict.items():
		tot3 += v
		for K,V in scoringdict1.items():
			if k == K:
				mut_intrascore += (v*V) #/len(mut_intradict)

	for k1,v1 in mut_interdict.items():
		tot4 += v
		for K1,V1 in scoringdict2.items():
			if k1 == K1:
				mut_interscore += (v1*V1) #/len(mut_interdict)



#	print('wt_interscore = ',wt_intrascore)
#	print('wt_intrascore = ',wt_interscore)
#	print('mut_intrascore = ',mut_intrascore)
#	print('mut_interscore = ',mut_interscore)

#	print(wt_intrascore+wt_interscore)

	DDG_pred = ((mut_intrascore)+(mut_interscore)) - ((wt_intrascore)+(wt_interscore))
#	print(mut)
	print(mut+'   '+'DDG_pred ='+' '+str(DDG_pred))



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
#		if 'Pr/PI' not in lines or 'AB/AG' not in lines:
		lines = lines.split(',')
		mutation = mut_converter(lines)
		mut_list.append(mutation)
		pdb_id = lines[0]
		ids = pdb_id.split('_')
		id = ids[0]
		if id not in id_list:
			id_list.append(id)

	for mut in mut_list:
		pdb_id = mut.split('_')[0]
		chain = mut.split(' ')[2]
		wt_res = mut.split(' ')[1]
		mut_position = mut.split(' ')[3]
		mut_res = mut.split(' ')[4]
		intercontacts = interface(pdb_id)
		filt_relacc_intercontacts = filt_relsurf_acc(intercontacts,pdb_id)
		intracont = intracontacts(pdb_id)
		prediction(matrix1,matrix2,filt_relacc_intercontacts,intracont,mut)
#		break

