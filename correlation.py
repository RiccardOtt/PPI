import sys
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from numpy import mean
from numpy import std
from numpy import cov
from scipy.stats.stats import pearsonr
from numpy.random import randn
from numpy.random import seed
from matplotlib import pyplot
import matplotlib.pyplot as plt
from numpy.random import randn
from numpy.random import seed


def correlation(ext,mine):

#	extl = ext[np.tril_indices(20)]
#	extu = ext[np.triu_indices(20)]
#	minel = mine[np.tril_indices(20)]
#	mineu = mine[np.triu_indices(20)]

#	corr = np.corrcoef(extl,minel)
#	print(corr)

	empty = np.zeros((20,20))

	A = np.tril(mine) + np.triu(mine.T,1)
	A = pd.DataFrame(A)
	print(A)

#	m = pd.DataFrame(m)
	A.to_csv('BASU.txt',sep=' ')



if __name__ == '__main__':
	skol_matrix = sys.argv[1]
	my_matrix = sys.argv[2]

	skolnic = np.loadtxt(skol_matrix)
	matrix = np.loadtxt(my_matrix)

#	skolnic = pd.read_csv(skol_matrix, delimiter = "\t")
#	matrix = pd.read_csv(my_matrix,delimiter = "\t")
#	matrix = np.array(matrix)

	correlation(skolnic,matrix)
