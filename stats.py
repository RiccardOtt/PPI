import sys
import numpy as np
import re
from numpy import random
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.stats import norm
import os
import filecmp



def divide_singl_to_mult_mut(skempi):
	ids_single_mut = []
	ids_mult_mut = []
	single_mut = []
	database = open(skempi,'r')

	for lines in database:
		lines = lines.split(';')
		if ',' in lines[1]:
			continue
		else:
			ids_single_mut.append(lines)

	return ids_single_mut,ids_mult_mut



def mut_location(db,ids):
	sup = 0
	cor = 0
	intt = 0
	sur = 0
	rim = 0
	database = open(db,'r')
	s = []
	c = []
	i = []
	su = []
	r = []

	for lines in database:
		lines = lines.rstrip().split(';')
		if lines[3] == 'SUP':
			sup += 1
			s.append(lines)
		elif lines[3] == 'COR':
			cor += 1
			c.append(lines)
		elif lines[3] == 'INT':
			intt += 1
			i.append(lines)
		elif lines[3] == 'SUR':
			sur += 1
			su.append(lines)
		else:
			rim += 1
			r.append(lines)

	labels = ['sup','cor','int','sur','rim']
	colors = ['b','c','r','g','m','y']
	sizes = [sup,cor,intt,sur,rim]

	fig1, ax1 = plt.subplots()
	ax1.pie(sizes, colors=colors, labels=labels, autopct='%1.1f%%', startangle=90)
	centre_circle = plt.Circle((0,0),0.70,fc='white')
	fig=plt.gcf()
	fig.gca().add_artist(centre_circle)
	ax1.axis('equal')
	plt.title('Skempi mutations localization')
#	plt.show()

	return s,c,i,su,r


def take_values(db,ids):
	ln = np.log
	R = float(0.001985)
	db = open(db)
	ids_singlemut_filt = []
	mut = []
	affinity_wt = []
	affinity_mut = []
	delta_G_wt = []
	delta_G_mut = []
	localization = []
	db_singlmut = []


	for lines in ids:
		if 'n' in lines[6] or 'n' in lines[8] or '>' in lines[6] or '>' in lines[8] or '<' in lines[6] or '<' in lines[8] or '~' in lines[6] or '~' in lines[8] or not lines[13]:
			continue
		else:
#			print(lines)
			db_singlmut.append(lines)
			ids_singlemut_filt.append(lines[0])
			mut.append(lines[2])
			affinity_wt.append(float(lines[8]))
			affinity_mut.append(float(lines[6]))
			delta_G_wt.append(float(re.sub('[^0-9]','',lines[13]))*R*ln(float(lines[8])))
			delta_G_mut.append(float(re.sub('[^0-9]','',lines[13]))*R*ln(float(lines[6])))
			localization.append(lines[3])

	for i in range(len(db_singlmut)):
		print(ids_singlemut_filt[i]+' '+mut[i]+' '+str(delta_G_mut[i]-delta_G_wt[i]))

	return affinity_wt,affinity_mut,delta_G_wt,delta_G_mut,set(ids_singlemut_filt),db_singlmut,localization



def stats(wt,mut,d_G_wt,d_G_mut,skem_singmut,local):
	c_wt = []
	c_mut = []

	WT = []
	MUT = []

	Delta_Delta_G = []

	stab_mut = []
	destab_mut = []

	val_wt = np.array(wt)
	val_mut = np.array(mut)
	val_dG_wt = np.array(d_G_wt)
	val_dG_mut = np.array(d_G_mut)

	mean_wt = np.mean(val_wt)
	mean_mut = np.mean(val_mut)
	sd_wt = np.std(val_wt)
	sd_mut = np.std(val_mut)

	mean_dG_wt = np.mean(val_dG_wt)
	mean_dG_mut = np.mean(val_dG_mut)
	sd_dG_wt = np.std(val_dG_wt)
	sd_dG_mut = np.std(val_dG_mut)

#	print('Mean aff wt =', mean_wt)
#	print('Mean aff mut =', mean_mut)
#	print('stand. aff dev. wt =', sd_wt)
#	print('stand. aff dev. mut =', sd_mut)

#	print('Mean deltaG wt =', mean_dG_wt)
#	print('Mean deltaG mut =', mean_dG_mut)
#	print('stand. dev. deltaG wt =', sd_dG_wt)
#	print('stand. dev. deltaG mut =', sd_dG_mut)


	for i in range(len(d_G_wt)):
		Delta_Delta_G.append((d_G_mut[i] - d_G_wt[i]))
#	print(Delta_Delta_G)

#	for i in range(len(Delta_Delta_G)):
#		print(local[i], Delta_Delta_G[i])


	for i in range(len(Delta_Delta_G)):
		if Delta_Delta_G[i] <= 0:
#			print(len(skem_singmut[i]))
			stab_mut.append(Delta_Delta_G[i])
		if Delta_Delta_G[i] > 0:
#			print(len(skem_singmut[i]))
			destab_mut.append(Delta_Delta_G[i])


	DDG = np.array(Delta_Delta_G)
	mean_ddG = np.mean(DDG)
	sd_DDG = np.std(DDG)
#	print('Mean DDG =', mean_ddG)
#	print('stand. dev. DDG =', sd_DDG)

	destab_mut = np.array(destab_mut).astype(np.float)
	stab_mut = np.array(stab_mut).astype(np.float)
#	print('mean destabilizing mutatios =', np.mean(destab_mut))
#	print('mean stabilizing mutations =', np.mean(stab_mut))
#	print('std destabilizing mutatios =', np.std(destab_mut))
#	print('std stabilizing mutatios =', np.std(destab_mut))


	np.set_printoptions(suppress=True)


######  PLOT THE DISTRIBUTION #########
#	fig = sns.kdeplot(destab_mut,shade=True,color='r',label="Destabilizing")
#	fig = sns.kdeplot(stab_mut, shade=True,color='b',label="Stabilizing")
#	fig = sns.kdeplot(DDG,shade=True,color='b')
#	plt.title('Normal Distributions DeltaDeltaG')
#	plt.legend(labels=['DeltaDeltaG'], loc='upper right')
#	plt.xlabel('Kcal/mol')
#	plt.show()

#	print(len(stab_mut)/len(DDG))



def take_chain(idschains,pdb_fl):

	pdb_lines = ''
	idschains = idschains.split('_')
	chains = idschains[1:]
	chains = chains[0]
	for lines in pdb_fl:
		if 'TER' not in lines:
			if lines[21] == chains:
				pdb_lines += lines

	return pdb_lines




def parse_dssp(dssp_file,id,mutsplit,loc):
	res_acc = []
	for line in dssp_file:
		line = line.rstrip()
		dssp_target = line[5:15]
		reschain_mut = mutsplit[0]
		mut_position = mutsplit[1]
		if reschain_mut[0] in dssp_target:
			if reschain_mut[1] in dssp_target:
				if mut_position in dssp_target:
					res_acc.append(str(id)+' '+reschain_mut+' '+mut_position+' '+line[35:38]+' '+loc)

	return res_acc




if __name__ == '__main__':
	db = sys.argv[1]
	singmut = []
	ids_single, ids_multiple = divide_singl_to_mult_mut(db)
#	s,c,i,su,r = mut_location(db,ids_single)
	aff_wt, aff_mut,del_G_wt,del_G_mut,pdb_ids_singlemut_filt,skem_singlmut,localization = take_values(db,ids_single)
	stats(aff_wt, aff_mut,del_G_wt,del_G_mut,skem_singlmut,localization)


#######MAKE DSSP##############
#	for i in skem_singlmut:
#		singmut.append(i[2])

#	set_id = set()       ####################################
#	set_id_chain = set() ###for not take chain multiple times
#	id_and_chain_mut = []
#	for lines in skem_singlmut:
#		idschains = lines[0].split('_')
#		id = idschains[0]
#		set_id.add(id)
#		chains = lines[0].split('_')[1:]
#		chain = ''.join(chains)
#		chain_split = list(chain)
#		for ch in chain_split:
#			set_id_chain.add(id+'_'+ch)

#	for i in set_id_chain:
#		g = i.split('_')[0]
#		pdb_file = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/'+g+'.pdb','r')
#		pdb_chain = take_chain(i,pdb_file)
#		pdbchain = open('/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/onlymut_chain/'+i+'.pdb','w')
#		pdbchain.write(pdb_chain)
#		pdbchain.close()

#	for el in set_id:
#		os.system('mkdssp -i' +' '+ '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/complex/' + el + '.pdb' +' '+'-o' +' '+ '/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/complex/' + el + '.dssp')
#	for el in set_id_chain:
#		os.system('mkdssp -i' +' '+ '/home/riccardo/Documents/Bioinf/Documents/Tesi/pdb/PDBs/onlymut_chain/' + el + '.pdb' +' '+'-o' +' '+ '/home/riccardo/Documents/Bioinf/Documents/Tesi/dssp_single_mut/onlymut_chain/' + el + '.dssp')



########SINGLE CHAIN REL SURF ACC######
#	surf_acc = []
#	id_chain_mut_toopen = []
#	for lines in skem_singlmut:
#		pdb_id = lines[0]
#		ids = pdb_id.split('_')
#		id = ids[0]
#		chains = ids[1:]
#		print(id)
#		chain = ''.join(chains)
#		chain_split = list(chain)
#		mut = lines[2]
#		mut_split = re.split('(\d+)',mut)
#		res_chain = mut_split[0]
#		chain = res_chain[1]
#		loc = lines[3]
#		try:
#			dssp_file = open('/home/riccardo/Documents/Documents/Tesi/dssp_single_mut/onlymut_chain/'+id+'_'+chain+'.dssp','r')
#		except:	continue
#		else:
#			acc = parse_dssp(dssp_file,id,mut_split,loc)
#			surf_acc.append(acc[0])

#	print(surf_acc)
#	for i in surf_acc:
#		print(i)

#	surfacc_array = np.array(surf_acc).astype(np.float)
#	meansurf = np.mean(surfacc_array)
#	stddev_surf = np.std(surfacc_array)
#	print('mean surface access single mutated chain =', meansurf)
#	fig = sns.kdeplot(surfacc_array,shade=True,color='r')
#	plt.show()



########COMPLEXES RELATIVE SURF ACC FOR MUTATED RES#########
#	relsurf_acc = []
#	for lines in skem_singlmut:
#		pdb_id = lines[0]
#		ids = pdb_id.split('_')
#		id = ids[0]
#		mut = lines[2]
#		mut_split = re.split('(\d+)',mut)
#		dssp_file = open('/home/riccardo/Documents/Documents/Tesi/dssp_single_mut/complex/'+id+'.dssp','r')
#		acc = parse_dssp(dssp_file,id,mut_split,loc)
#		relsurf_acc.append(acc[0])
#		loc = lines[3]

#	print(relsurf_acc)
#	for i in relsurf_acc:
#		print(i)

#	surfacc_array = np.array(surf_acc).astype(np.float)
#	meansurf = np.mean(surfacc_array)
#	stddev_surf = np.std(surfacc_array)
#	print('mean surface access complex =', meansurf)
#	fig = sns.kdeplot(surfacc_array,shade=True,color='r')
#	plt.show()


#######COMPARING REL SURF ACC COMPLEX AND SINGLE CHAIN########
####### bash on the two files produces: rel_sol...complex.txt, rel_sol....singlchain.txt##########
