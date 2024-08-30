![](amapec_logo.svg)

## About
AMAPEC is a predictor of antimicrobial activity for fungal secreted proteins, that aims to assist researchers in the characterization of new effectors. By offering unprecedented functional insights into fungal secretomes (generally sparsely functionally annotated), it may aid in biological interpretations during genomic, transcriptomic and proteomic analyses.

Using a (predicted) protein structure as input, AMAPEC returns:
- a mean confidence score for the input predicted protein structure (pLDDT⁠), with the rationale that a low-confidence structure may obtain a predicted antimicrobial activity that is not biologically meaningful
- a classification as ‘Antimicrobial’ or ‘Non-antimicrobial’
- a probability score for its antimicrobial activity, that ranges between 0 (no antimicrobial activity) and 1 (highly likely to be antimicrobial).

#### Citation

Mesny, F. & Thomma, B. P. (2024). AMAPEC: accurate antimicrobial activity prediction for fungal effector proteins. *BioRxiv*, 2024-01.
doi: [10.1101/2024.01.04.574150](https://www.biorxiv.org/content/10.1101/2024.01.04.574150)

## Use AMAPEC online

[Click here](https://colab.research.google.com/github/fantin-mesny/amapec/blob/main/googleColab/AMAPEC.ipynb) to try AMAPEC online and compute antimicrobial activity predictions using Google Colab.

## Installation

### Main installation

You can simply download this repository or clone it to your current directory:
```
git clone https://github.com/fantin-mesny/amapec
chmod +x amapec/amapec
```
### Dependencies

We recommend to install AMAPEC dependencies with `conda`, using the file environment.yml:
```
conda env create -f environment.yml
```
If the `conda` installation fails to resolve dependencies, we recommend trying to install AMAPEC with `mamba`.

(facultative) If you wish to re-train the prediction model yourself, additional installation of software [mercat2](https://github.com/raw-lab/mercat2) is needed.
Please follow the developer's installation guidelines and install the software in a different mamba environment named `mercat2`:
```
conda create -n mercat2 -c conda-forge -c bioconda mercat2
```
### Check installation

Run a test prediction on the provided test set:
```
conda activate amapec
amapec -i testSet 
```
See the output file in directory testSet_AMprediction to check if the prediction worked.
The log.txt file in the output directory should inform you about the origins of eventual problems.

## Running AMAPEC

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

We provide in the repository the full training dataset as a zipped archive as well as scripts to retrain the SVM classifier of AMAPEC.
Information regarding the composition of the training dataset can be found in and in the AMAPEC preprint (see citation above) and in the associated supplementary tables.

```
conda activate amapec
bash trainModel.sh
```
You can edit the script header to modify inputs (positive and negative training sets) and parameters (e.g. number of threads).  
Output from this training pipeline will be written in folder "retrained_model" in the amapec directory.

You can then run amapec with the retrained model in the following way:
```
amapec -i /path/to/input_dir -o /path/to/output_dir -m retrained_model/model -t 40
```

## Recommendations

We highly recommend to use AMAPEC on protein structures that were predicted from mature sequences (after signal peptide removal). 
Structure prediction with a signal peptide can considerably affect protein geometries, and subsequently, impact antimicrobial activity prediction.
For this reason, using effector structures downloaded from the AlphaFold database is not recommended.

Since AlphaFold is computationally demanding and may be difficult to run on a large number of proteins, we recommend the use of [ESM-Fold](https://github.com/facebookresearch/esm) or alternatively, of [ColabFold](https://github.com/sokrypton/ColabFold). 

## References

AMAPEC reimplements portions of code from the following repositories:
- [harmslab/pdbtools](https://github.com/harmslab/pdbtools) (GPL-3.0 license)
- [ugSUBMARINE/strucural-properties](https://github.com/ugSUBMARINE/structural-properties) (MIT license)
- [sarisabban/Rg](https://github.com/sarisabban/Rg) (MIT license)
- [SBRG/ssbio](https://github.com/SBRG/ssbio) (MIT license)

The full list of references can be found in the AMAPEC preprint (reference above), in the **References** section and in Supplementary Table 4.

## Contact

fmesny1 \[at\] uni-koeln.de
