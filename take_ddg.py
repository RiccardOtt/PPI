import sys
import numpy as np
import re


def take_values(db):
	ln = np.log
	R = float(0.001985)
	ids_singlemut_filt = []
	affinity_wt = []
	affinity_mut = []
	delta_G_wt = []
	delta_G_mut = []
	localization = []
	db_singlmut = []


	for lines in db:
		if 'n' in lines[6] or 'n' in lines[8] or '>' in lines[6] or '>' in lines[8] or '<' in lines[6] or '<' in lines[8] or '~' in lines[6] or '~' in lines[8] or not lines[13]:
			continue
		else:
			ids_singlemut_filt.append(lines[0])
			affinity_wt.append(float(lines[8]))
#			affinity_mut.append(float(lines[6]))
#			print(affinity_wt)
#			delta_G_wt.append(float(re.sub('[^0-9]','',lines[13]))*R*ln(float(lines[8])))
#			delta_G_mut.append(float(re.sub('[^0-9]','',lines[13]))*R*ln(float(lines[6])))
#			localization.append(lines[3])

#	return affinity_wt,affinity_mut,delta_G_wt,delta_G_mut,set(ids_singlemut_filt),db_singlmut,localization



if __name__ == '__main__':
	skempi = sys.argv[1]

	with open(skempi) as db:
		take_values(db)
