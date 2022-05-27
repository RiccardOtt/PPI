import sys
import re
import Bio.PDB
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB import *




def interface(id):

	amino = ['GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP']
	contacts = []
	parser = PDBParser()

	structure = parser.get_structure('PHA-L', '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/'+id+'.pdb')

#	for model in structure:
#		for chain in model:
#			for residue in chain:
#				print(residue)

	model = structure[0]
	for chain1 in model:
		for chain2 in model:
			if chain1 != chain2:
				for residue1 in chain1:
					tags1 = residue1.id
#					print(residue1)
					for residue2 in chain2:
						tags2 = residue2.id
#						if residue1 != residue2:
						if tags1[0] != " " or tags2[0] != " ":
					  	  ###	The residue is a heteroatom
							pass
						else:
#							try:
							atoms1 = residue1.get_atoms()
							atoms2 = residue2.get_atoms()
							for atom1 in atoms1:
								for atom2 in atoms2:
									distance = atom1 - atom2
#									print(distance)
#							distance = residue1['CA'] - residue2['CA']
#									except KeyError:
        	       				## no CA atom, e.g. for H_NAG
#									continue
									if distance < int(6):
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


	for i in cont_filt:
		print(i)










if __name__ == '__main__':
	skempi_single = sys.argv[1]
	skem = open(skempi_single)
	id_list = []


	for lines in skem:
		if 'Pr/PI' in lines or 'AB/AG' in lines:
			lines = lines.split(',')
			pdb_id = lines[0]
			ids = pdb_id.split('_')
			id = ids[0]
			if id not in id_list:
				id_list.append(id)


	for pdb_id in id_list:
		intercontacts = interface(pdb_id)
		filt_relacc_intercontacts = filt_relsurf_acc(intercontacts,pdb_id)
#		break
