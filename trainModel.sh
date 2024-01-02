#!/bin/bash

scriptDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
positiveSet=$scriptDir/positiveSet
negativeSet=$scriptDir/negativeSet
workingDir=$scriptDir/retrained_model
nproc=40
mkdir -p $workingDir

# Define function checking intermediary output files
function check_output () {
  local FILE=$1
  if [ -f $FILE ]; then
    echo "File" $FILE "has been written." >> $workingDir/log.txt 2>&1
  else
    echo "ERROR. File" $FILE "cannot be found. Please check the errors in log file" $workingDir/log.txt
    exit 1
  fi
}



############ Calculate Sequence Properties ##################
echo "Calculating sequence properties..."
cat $positiveSet/sequences/* $negativeSet/sequences/* > $workingDir/sequences.fasta
Rscript $scriptDir/scripts/calculateSequenceProperties.R $workingDir/sequences.fasta $workingDir/seqProperties.csv >> $workingDir/log.txt 2>&1
check_output $workingDir/seqProperties.csv

############ Get Kmer content of sequences ##################
echo "Translating sequences into reduced amino-acid alphabet..."
python $scriptDir/scripts/transformSeqs.py -i $workingDir/sequences.fasta -o $workingDir/transformedSeqs.fa >> $workingDir/log.txt 2>&1
check_output $workingDir/transformedSeqs.fa

echo "Profiling K-mer content..."
mkdir -p $workingDir/mercat
conda run -n mercat2 mercat2.py -i $workingDir/transformedSeqs.fa -k 3 -o $workingDir/mercat/mercat_output_k3 >> $workingDir/log.txt 2>&1 #### install mercat in env
conda run -n mercat2 mercat2.py -i $workingDir/transformedSeqs.fa -k 4 -o $workingDir/mercat/mercat_output_k4 >> $workingDir/log.txt 2>&1
conda run -n mercat2 mercat2.py -i $workingDir/transformedSeqs.fa -k 5 -o $workingDir/mercat/mercat_output_k5 >> $workingDir/log.txt 2>&1
conda run -n mercat2 mercat2.py -i $workingDir/transformedSeqs.fa -k 6 -o $workingDir/mercat/mercat_output_k6 >> $workingDir/log.txt 2>&1
python $scriptDir/scripts/mercatOutputParser.py -s $workingDir/transformedSeqs.fa -m $workingDir/mercat -o $workingDir/kmerContent.csv >> $workingDir/log.txt 2>&1
check_output $workingDir/kmerContent.csv

########### Calculate Structure Properties ##################
echo "Calculating structure properties..."
python $scriptDir/scripts/calculateStructureProperties.py -nproc $nproc -i $positiveSet/structures/,$negativeSet/structures/ -o $workingDir >> $workingDir/log.txt 2>&1 
check_output $workingDir/strucProperties.csv

###################### Build model ##########################
echo "Training model and performing cross-validation..."
mkdir -p $workingDir/model
python $scriptDir/scripts/trainSVM.py --svmScale -o $workingDir/model -seq $workingDir/seqProperties.csv -struc $workingDir/strucProperties.csv -kmer $workingDir/kmerContent.csv -posPrefix a,r,o -negPrefix n,q,m --weighted -posMeta $positiveSet/metadata.csv -negMeta $negativeSet/metadata.csv >> $workingDir/log.txt 2>&1
check_output $workingDir/model/classifier.sav

