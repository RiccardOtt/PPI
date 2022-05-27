import sys
import pandas as pd
import numpy as np
import re


def get_matrix(aaindexfile,matrixname):
	canPrintLines = False
	for line in aaindexfile:
		if matrixname in line:
			canPrintLines = True # We have found an @ so we can start printing lines
		elif '//' in line:
			canPrintLines = False # We have found a + so we don't want to print anymore
		if canPrintLines:
			line = line.rstrip()
			line = line[2::]
			if re.match(r'\W', line):
#				arr = np.array(line)
				print(line)



if __name__ == '__main__':
	aaindex = sys.argv[1]
	getmatrix = sys.argv[2]
	aa = open(aaindex)
	get_matrix(aa,getmatrix)
