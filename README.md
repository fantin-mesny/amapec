# amapec v1.0

## About
amapec is...

## Installation

### Dependencies

We recommend installing amapec dependencies with `mamba`, using the file environment.yml:
```
mamba env create -f environment.yml
```
Installation with `conda` should also be possible but might fail to resolve dependencies.

If you wish to re-train the prediction model yourself, additional installation of software [mercat2](https://github.com/raw-lab/mercat2) is needed.
Please follow the developer's installation guidelines and install the software in a different mamba environment named `mercat2`:
```
mamba create -n mercat2 -c conda-forge -c bioconda mercat2
```

### Main installation

You can simply download this repository or clone it to your current directory:
```
git clone https://github.com/fantin-mesny/amapec
```
### Check installation

Run a test prediction on the provided test set:
```
conda activate amapec
amapec -i testSet 
```
See the output file in directory testSet_AMprediction to check if the prediction worked.
The log.txt file in the output directory should inform you where the problem originates.

## Running amapec

### Simple prediction

To predict antimicrobial activity of a set of proteins which PDB structures are in directory "input", using 40 CPU threads:
```
conda activate amapec
amapec -i input/ -o amapec_prediction -t 40 
```
In the "amapec_prediction" directory, you will find a "prediction.csv" file containing the prediction results.

### Command line arguments

```
amapec [-i <directory>] [-o <directory>] [-t <num>] [-h]
  -i: directory containing protein structures in PDB format (required)
  -t: number of threads/cores to use (default: 4)
  -o: output directory (default: <input_directory>_AMprediction)
  -d: with this option specified, temporary files will not be deleted (for debugging purposes).
  -h: display help
```
### Re-training the model

We provide in the repository the full training dataset as well as scripts to retrain the SVM classifier amapec is based on.
```
conda activate amapec
bash trainModel.sh
```
Output from this training pipeline will be written in folder "retrained_model" in the amapec directory.

## Contact

fmesny1 \[at\] uni-koeln.de
