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
			for residue1 in chain1:
				tags1 = residue1.id
#				print(residue1)
				for residue2 in chain2:
					tags2 = residue2.id
#					if residue1 != residue2:
					if chain1 != chain2:
						if tags1[0] != " " or tags2[0] != " ":
						    ### The residue is a heteroatom
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
				print(lines[0:7],lines[9:16])
		except:
			pass






if __name__ == '__main__':

	skempi_singlmut = sys.argv[1]

	set_id_chain = set()
	set_id = set()
	skempi_singl = open(skempi_singlmut)

	for lines in skempi_singl:
		lines = lines.split(',')
		pdb_id = lines[0]
		ids = pdb_id.split('_')
		id = ids[0]
#		chains = ids[1:]
#		chain = ''.join(chains)
#		chain_split = list(chain)
#		for ch in chain_split:
#			set_id_chain.add(id+'_'+ch)
		set_id.add(id)


	for i in set_id:
#	for i in set_id_chain:
#		g = i.split('_')[0]
#		try:
#			dssp_file_complex = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/complex/'+g+'.dssp','r')
#			dssp_file_singlchain = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/onlymut_chain/'+i+'.dssp','r')
#		except: continue
#		else:
		cont = interface(i)
		contacts_filtered = filt_relsurf_acc(cont,i)
#		break
