import os
import sys
import argparse
import pandas as pd
from sklearn import svm
import joblib
from Bio import SeqIO

def get_params(argv):
	parser = argparse.ArgumentParser(description='Use pre-computed model to predict antimicrobial activity from a CSV input with protein properties')
	parser.add_argument('-m', '--m', help="Classifier file (.sav)", required=True)
	parser.add_argument('-pe', '--pe', help="Probability estimator file (.sav)", required=True)
	parser.add_argument('-struc_prop','--struc_prop',help='Input CSV with structural properties of proteins to classify',required=True)
	parser.add_argument('-seq_prop','--seq_prop',help='Input CSV with sequence properties of proteins to classify',required=True)
	parser.add_argument('-seq','--seq',help='Fasta file with transformed sequences',required=True)
	parser.add_argument('-o', '--o', help="Output directory", required=True)
	a = parser.parse_args()
	return a


if __name__ == '__main__':
	a = get_params(sys.argv[1:])

	## load the model from disk
	classifier = joblib.load(a.m)
	prob_estimator=joblib.load(a.pe)
	model_features=classifier.feature_names_in_
	
	## load protein physicochemical properties
	seq=pd.read_csv(a.seq_prop).set_index('Name').drop(columns='Unnamed: 0').sort_index() #load sequence properties
	struc=pd.read_csv(a.struc_prop).set_index('Unnamed: 0').sort_index() # load structure properties
	struc=struc.astype(float)
	with open(a.seq,'r') as fa:
		fasta=SeqIO.to_dict(SeqIO.parse(fa,'fasta')) # load transform sequences
	
	## counting k-mers
	kmerToCount=['0051','0310','1013','1115','6100','444'] #k for k in model_features if k not in seq and k not in struc]
	kmer_counts=pd.DataFrame(index=list(fasta.keys()), columns=list(kmerToCount))
	for s in kmer_counts.index:
		for kmer in kmer_counts.columns:
			kmer_counts.loc[s,kmer]=fasta[s].seq.count(kmer)
	kmer_counts=kmer_counts.sort_index()

	## merge protein data into one dataframe
	df=pd.concat([seq,struc,kmer_counts],axis=1).dropna()
	df=df[[c for c in model_features]]

	## transform test dataset to be in the same scale as the training one
	scale=pd.read_csv(a.m+'.scale').set_index('Unnamed: 0')
	for column in df.columns: 
		df[column]=(df[column] - scale.loc[column,'Mean']) / scale.loc[column,'Standard deviation']

	## run prediction
	pred=pd.DataFrame(prob_estimator.predict_proba(df),columns=['Probability of class '+str(int(cla)) for cla in prob_estimator.classes_], index=df.index)
	predi=pd.DataFrame(classifier.predict(df),columns=['Prediction'], index=df.index)
	
	## prepare output dataframe
	pLDDTs=pd.read_csv(a.struc_prop.replace('strucProperties.csv','pLDDTs.csv')).set_index('Unnamed: 0').sort_index()
	pred=pd.concat([pLDDTs,pred,predi],axis=1)
	pred['Prediction']=pred['Prediction'].map({0:'Non-antimicrobial',1:'Antimicrobial'})
	pred=pred.drop(columns=['Probability of class 0']).rename(columns={'Probability of class 1':'Probability of antimicrobial activity'})
	pred.index.name='Protein ID'
	
	if os.path.isfile(a.o+'/renamedFiles.tsv'):
		rename=pd.read_csv(a.o+'/renamedFiles.tsv',sep='\t',header=None).set_index(1)[0].to_dict()
		pred.index=pred.index.map(rename)
	
	pred.to_csv(a.o+'/prediction.csv')


