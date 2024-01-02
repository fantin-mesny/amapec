import sys
from Bio import SeqIO
import argparse
import os
import subprocess
import multiprocessing
from Bio import BiopythonParserWarning
import warnings

def get_params(argv):
	parser = argparse.ArgumentParser(description='Create sequence fasta file from structures')
	parser.add_argument('-i', '--i', help="Directory including structures'", required=True)
	parser.add_argument('-o', '--o', help="Directory to output sequences.fasta in'", required=True)
	parser.add_argument('-nproc', '--nproc', help="Nproc", default=50)
	a = parser.parse_args()
	return a

def fastaFromPDB(pdb,a,return_list):
	with open(a.i+'/'+pdb, 'r') as pdb_file:
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', BiopythonParserWarning)
			record=list(SeqIO.parse(pdb_file, 'pdb-atom'))
	if len(record)==1:
		record[0].name=pdb.replace('.pdb','')
		record[0].id=pdb.replace('.pdb','')
		record[0].description=pdb.replace('.pdb','')
	else:
		print('ERROR: multiple structure found in pdb file', pdb)
	return_list.append(record[0])

def sequenceContainsX(seqRecord):
	if 'X' in str(seqRecord.seq):
		return True
	else:
		return False
		
def filterOutSequencesWithX(return_list):
	toIgnore=[]
	return_list_filtered=[]
	for seqRecord in return_list:
		if sequenceContainsX(seqRecord):
			toIgnore.append(seqRecord.id)
		else:
			return_list_filtered.append(seqRecord)
	
	if len(toIgnore)>0:
		with open(a.o+'/ignored.txt','w+') as output:
			for seqId in toIgnore:
				output.write(seqId+'\n')			
	return return_list_filtered


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	pool = multiprocessing.Pool(int(a.nproc))
	manager = multiprocessing.Manager()
	return_list = manager.list()

	for pdb in os.listdir(a.i):
		if pdb.endswith('.pdb'):
			pool.apply_async(fastaFromPDB, args=(pdb,a,return_list))
	pool.close()
	pool.join()

	return_list=filterOutSequencesWithX(return_list)
	
	with open(a.o+'/sequences.fasta','w+') as outp:	
		SeqIO.write(return_list,outp,'fasta')
