import sys
import os


def stat(mut,db):
	for lines in db:
		lines = lines.split(';')
		if mut == lines[1]:
			print(lines)









if __name__ == '__main__':
	skempi = sys.argv[1]
	mutation_repeated = sys.argv[2]
	mut = open(mutation_repeated)
	db = open(skempi)
	for i in mut:
		stat(i,db)
