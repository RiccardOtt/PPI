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


def correlation(skol,mine):

	upper_skol = skol[np.triu_indices(20)]
	upper_mine = mine[np.triu_indices(20)]
#	skol = skol.drop(['Unnamed: 20'], axis = 1)
#	skol = np.array(skol)
#	mine = np.array(mine)
#	print(skol)
	corr = np.corrcoef(upper_skol, upper_mine)
	print(corr)

#	s = pd.DataFrame(skol)
#	m = pd.DataFrame(mine)
#	print(skol)
#	print(mine)
#	print(s.corrwith(other=m))




if __name__ == '__main__':
	skol_matrix = sys.argv[1]
	my_matrix = sys.argv[2]
	skolnic = np.loadtxt(skol_matrix)
	matrix = np.loadtxt(my_matrix)
#	skolnic = pd.read_csv(skol_matrix, delimiter = "\t")
#	matrix = pd.read_csv(my_matrix,delimiter = "\t")
	correlation(skolnic,matrix)
