from Bio import SeqIO
import argparse
import sys
import os

def get_params(argv):
	parser = argparse.ArgumentParser(description='Translate sequences from a file into a selected reduced aminoacid alphabet')
	parser.add_argument('-i', '--i', help="Fasta file containing all sequences'", required=True)
	parser.add_argument('-o', '--o', help="Output fasta file", required=True)
	parser.add_argument('-a', '--a', help="Alphabet to be used'", default='6')
	a = parser.parse_args()
	return a

alphabets={
	'1':'AG-C-DEKNPQRST-FILMVWY-H', # Accuracy of Sequence Alignment and Fold Assessment Using Reduced Amino Acid Alphabets (2006) #10
	'2':'ED-QST-NH-YPRK-LMI-VAF-GW-C', # Deep Learning Improves Antimicrobial Peptide Recognition (2018) #48
	'3':'FILVWY-ACGMP-DEHNR-KQST', # Reduced amino acid alphabet is sufficient to accurately recognize intrinsically disordered protein (2004) #18
	'4':'DE-KR-NQST-GP-HWYF-ACMLIV', # Universally Conserved Positions in Protein Folds:Reading Evolutionary Signals about Stability, Folding Kinetics and Function (1999) #37
	'5':'WFYLIVM-DE-KR-ACGHNPQST', # Predicting Disordered Regions in Proteins Based on Decision Trees of Reduced Amino Acid Composition (2006) #62
	'6':'GA-STC-VLIMP-FYW-DE-NQ-HKR' # small, nucleophilic,hydrophobic,aromatic,acidic, amide, basic https://international.neb.com/tools-and-resources/usage-guidelines/amino-acid-structures
}
# for the first 5, see Supplementary Data of https://doi.org/10.1016/j.csbj.2022.07.001

def parseAlph(alph):
	alph=alph.split('-')
	dic={}
	n=0
	for clust in alph:
		for letter in clust:
			dic[letter]=str(alph.index(clust))
	return dic

def transformSeq(seq,alphabet):
	alph=parseAlph(alphabet)
	for aa in seq:
		seq=seq.replace(aa,alph[aa])
	return seq
	
if __name__ == '__main__':
	a = get_params(sys.argv[1:])

	fa=list(SeqIO.parse(a.i,'fasta'))
	with open(a.o,'w+') as outp:
		for seq in fa:
			outp.write('>'+str(seq.name)+'\n'+transformSeq(str(seq.seq),alphabets[a.a])+'\n')

