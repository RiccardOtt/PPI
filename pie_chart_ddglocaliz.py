import sys
import matplotlib.pyplot as plt
import pandas as pd



def mut_localization(ddg_loc):
	sup_pos = 0
	sup_neg = 0
	cor_pos = 0
	cor_neg = 0
	int_pos = 0
	int_neg = 0
	sur_pos = 0
	sur_neg = 0
	rim_pos = 0
	rim_neg = 0




	for lines in ddg_loc:
		lines = lines.rstrip().split(' ')
		if lines[0] == 'SUP' and float(lines[1]) > 0:
			sup_pos += 1
		elif lines[0] == 'SUP' and float(lines[1]) <= 0:
			sup_neg += 1
		if lines[0] == 'COR' and float(lines[1]) > 0:
			cor_pos += 1
		elif lines[0] == 'COR' and float(lines[1]) <= 0:
			cor_neg += 1
		if lines[0] == 'INT' and float(lines[1]) > 0:
			int_pos += 1
		elif lines[0] == 'INT' and float(lines[1]) <= 0:
			int_neg += 1
		if lines[0] == 'SUR' and float(lines[1]) > 0:
			sur_pos += 1
		elif lines[0] == 'SUR' and float(lines[1]) <= 0:
			sur_neg += 1
		if lines[0] == 'RIM' and float(lines[1]) > 0:
			rim_pos += 1
		elif lines[0] == 'RIM' and float(lines[1]) <= 0:
			rim_neg += 1




	labels = ['sup stab','cor stab','int stab','sur stab','rim stab']
	colors = ['b','c','r','g','m','y']
	sizes = [sup_pos,cor_pos,int_pos,sur_pos,rim_pos]

	fig1, ax1 = plt.subplots()
	ax1.pie(sizes, colors=colors, labels=labels, autopct='%1.1f%%', startangle=90)
	centre_circle = plt.Circle((0,0),0.70,fc='white')
	fig=plt.gcf()
	fig.gca().add_artist(centre_circle)
	ax1.axis('equal')
	plt.title('Skempi stabilizing single mutations localization')
	plt.show()


	labels = ['sup destab','cor destab','int destab','sur destab','rim_destab']
	colors = ['b','c','r','g','m','y']
	sizes = [sup_neg,cor_neg,int_neg,sur_neg,rim_neg]

	fig1, ax1 = plt.subplots()
	ax1.pie(sizes, colors=colors, labels=labels, autopct='%1.1f%%', startangle=90)
	centre_circle = plt.Circle((0,0),0.70,fc='white')
	fig=plt.gcf()
	fig.gca().add_artist(centre_circle)
	ax1.axis('equal')
	plt.title('Skempi destabilizing single mutations localization')
	plt.show()




if __name__ == '__main__':
	ddglocalization = sys.argv[1]
	localiz = open(ddglocalization)
	mut_localization(localiz)
