
```
                                                _   ___  
   __ _ _ __ ___   __ _ _ __   ___  ___  __   _/ | / _ \ 
  / _` | '_ ` _ \ / _` | '_ \ / _ \/ __| \ \ / / || | | |
 | (_| | | | | | | (_| | |_) |  __/ (__   \ V /| || |_| |
  \__,_|_| |_| |_|\__,_| .__/ \___|\___|   \_/ |_(_)___/ 
                       |_|                               
                      
 AntiMicrobial Activity Prediction for Effector Candidates
   
```

## Installation

We recommend installing amapec with mamba, using the file environment.yml:
```
mamba env create -f environment.yml
```
Installation with conda is also possible but might fail to resolve dependencies.

If you wish to re-train the prediction model yourself, additional installation of software [mercat2](https://github.com/raw-lab/mercat2) is needed.
Please follow the developer's installation guidelines and install mercat2 in a different mamba environment:
```
mamba create -n mercat2 -c conda-forge -c bioconda mercat2
conda activate mercat2
```
