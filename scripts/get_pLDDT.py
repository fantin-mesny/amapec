import pandas as pd
import argparse
import sys
import os

def get_params(argv):
	parser = argparse.ArgumentParser(description='Calculate pLDDT from structures')
	parser.add_argument('-i', '--i', help="Comma-separated directories including structures'", required=True)
	parser.add_argument('-o', '--o', help="Output/working directory'", required=True)
	a = parser.parse_args()
	return a



def get_mean_pLDDT(File): #optimized for our AF2 structures
	with open(File,'r') as inp:
		lines=inp.readlines()
	lines_list=[l.split(' ') for l in lines]
	lines_list_clean=[[e for e in l if e!=''] for l in lines_list if l[0]=='ATOM']
	lines_list_filtered=[(str(l[4].replace('A','A,')+','+l[5]).replace(',,',',').split(',')[1],l[-3]) for l in lines_list_clean]
	vals=pd.DataFrame(lines_list_filtered)
	vals=vals.drop_duplicates()
	return vals[1].astype(float).mean()
	
	
if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	df=pd.DataFrame(columns=['pLDDT'])
	Dir=a.i+'/'
	
	if 'ignored.txt' in os.listdir(a.o):
		with open(a.o+'/ignored.txt') as inp:
			inpLines=inp.readlines()
		toIgnore=[ID.replace('\n','.pdb') for ID in inpLines]
	else:
		toIgnore=[]

	
	for struc in [st for st in os.listdir(Dir) if st.endswith('.pdb')]:
		if struc not in toIgnore:
			pLDDT=get_mean_pLDDT(Dir+struc)
			df.loc[struc.replace('.pdb',''),'pLDDT']=pLDDT
	df.sort_index().to_csv(a.o+'/pLDDTs.csv')
