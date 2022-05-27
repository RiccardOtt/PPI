import sys
from rosetta.protocols.scoring import Interface
from rosetta import *
from pyrosetta import *
#from pyrosetta.toolbox import cleanATOM
init()


def pyros_function(skempi_id):

	# 1. create a pose from the desired PDB file
	pose = pose_from_pdb('/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/1A22.pdb')
	dock_jump = 1
	interface = Interface(dock_jump)
	interface.distance(0.00000001)
	interface.calculate(pose)

	



if __name__ == '__main__':
	skempi_singlmut = sys.argv[1]
	set = set()
	skempi = open(skempi_singlmut)
	for lines in skempi:
		lines = lines.rstrip().split('_')
		id = lines[0]
		set.add(id)
	for ids in set:
		pyros_function(ids)
		break
