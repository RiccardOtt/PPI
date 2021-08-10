import sys
import os
from Bio.SeqUtils import seq3



def get_surface_res(dssp):

	relacc = {'G':104.0, 'A':129.0, 'S':155.0, 'C':167.0, 'P':159.0, 'T':172.0, 'V':174.0, 'N':195.0, 'D':193.0, 'L':201.0, 'I':197.0, 'M':224.0, 'E':223.0, 'Q':225.0, 'H':224.0, 'F':240.0, 'K':236.0, 'Y':263.0, 'W':285.0, 'R':274.0}
	for i, line in enumerate(dssp):
		if i > 27:
			line = line.rstrip()
			for k,v in relacc.items():
				if line[11] == k:
					rela = int(line[35:38])/int(v)*100
					if int(line.rstrip()[35:38]) > int(50):
						print(seq3(line[13]))







if __name__ == '__main__':
	for filename in os.listdir('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/only_chain/'):
		with open(os.path.join('/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/only_chain/', filename)) as f:
			get_surface_res(f)
