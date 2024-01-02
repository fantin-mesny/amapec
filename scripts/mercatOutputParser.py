import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO

def get_params(argv):
	parser = argparse.ArgumentParser(description='Parse mercat outputs, count detected kmer in sequences')
	parser.add_argument('-s', '--s', help="Fasta file containing all sequences", required=True)
	parser.add_argument('-m', '--m', help="Directory in which the mercat_output_* folders can be found", required=True)
	parser.add_argument('-o', '--o', help="Output CSV", required=True)
	a = parser.parse_args()
	return a


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	
	dfs_k=[]
	for folder in [folder for folder in os.listdir(a.m) if 'mercat_output_' in folder]:
		df_k=pd.read_csv(a.m+'/'+folder+'/combined_Nucleotide.tsv',sep='\t', dtype={'kmer': str})
		dfs_k.append(df_k)
	df_k=pd.concat(dfs_k).set_index('kmer')
	
	with open(a.s,'r') as fa:
		fasta=SeqIO.to_dict(SeqIO.parse(fa,'fasta'))
		
	df_o=pd.DataFrame(index=list(fasta.keys()), columns=list(df_k.index))
	for seq in df_o.index:
		for kmer in df_o.columns:
			df_o.loc[seq,kmer]=fasta[seq].seq.count(kmer)
	df_o.to_csv(a.o)
	print(df_o.sum(axis=0))
	
		
	
	
