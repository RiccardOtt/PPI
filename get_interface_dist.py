import sys
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)




def interface(id):

	contacts = []
#	set_int = set()
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
				for residue2 in chain2:
					if residue1 != residue2:
						if chain1 != chain2:
							try:
								distance = residue1['CA'] - residue2['CA']
							except KeyError:
              	       					## no CA atom, e.g. for H_NAG
								continue
							if distance < 4.5:
								contact = residue1,chain1,residue2,chain2,distance
#								set_int.add(str(residue1)+' '+str(residue2))
								contacts.append(contact)

	for i in contacts:
		print(i)

	return contacts




def get_residues(skempi):

	parser = PDBParser()

	structure = parser.get_structure('PHA-L', '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/'+id+'.pdb')

	for model in structure:
		for chain in model:
			for residue in chain:
				print(residue)






if __name__ == '__main__':

	skempi_singlmut = sys.argv[1]
	set_id = set()
	skempi = open(skempi_singlmut)
	for lines in skempi:
		lines = lines.rstrip().split('_')
		id = lines[0]
		set_id.add(id)
	for ids in set_id:
		cont = interface(ids)
#		get_residues(ids)
