#!/bin/bash

## Install miniconda
wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./miniconda.sh -b -f -p /usr/local
rm miniconda.sh

## Install mamba
conda config --add channels conda-forge
conda install -y mamba
mamba update -qy --all
mamba clean -qafy

## Install AMAPEC dependencies
mamba env create -f amapec/environment.yml
