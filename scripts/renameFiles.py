import re
import shutil
import argparse
import os
import sys

def get_params(argv):
	parser = argparse.ArgumentParser(description='Rename structure files')
	parser.add_argument('-i', '--i', help="Directory including structures'", required=True)
	parser.add_argument('-o', '--o', help="Directory to output renamed structures'", required=True)
	parser.add_argument('-nproc', '--nproc', help="Nproc", default=50)
	a = parser.parse_args()
	return a

def stringIncludesSpecialCharac(string):
	pattern = r"^[A-Za-z0-9._-]*$"
	return not re.match(pattern, string)


if __name__ == '__main__':
	a = get_params(sys.argv[1:])

	renamed={}
	includesSpec=0
	for f in os.listdir(a.i):
		if f.endswith('.pdb'):
			if stringIncludesSpecialCharac(f):
				includesSpec+=1
				#if ' ' in f:
				#	f2=f.split(' ')[0]+'.pdb'
				#else:
				f2=f
				
				for letter in f:
					if stringIncludesSpecialCharac(letter):
						f2=f2.replace(letter,'_')
				renamed[f]=f2
			else:
				renamed[f]=f
				
	if includesSpec>0:
		with open(a.o+'/renamedFiles.tsv','w+') as tsv:
			for p in renamed:
				tsv.write(p.replace('.pdb','')+'\t'+renamed[p].replace('.pdb','')+'\n')
				shutil.copyfile(a.i+'/'+p, a.o+'/structures/'+renamed[p])

	print(includesSpec)
