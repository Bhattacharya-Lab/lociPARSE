## lociPARSE: a locality-aware invariant point attention model for scoring RNA 3D structures

by Sumit Tarafder and Debswapna Bhattacharya

[[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1)] [[pdf](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1.full.pdf)]

Codebase for our <ins>loc</ins>ality-aware <ins>i</ins>nvariant <ins>P</ins>oint <ins>A</ins>ttention-based <ins>R</ins>NA <ins>S</ins>cor<ins>E</ins>r (lociPARSE).

![alt text](https://github.com/Bhattacharya-Lab/lociPARSE/blob/main/lociPARSE.png?raw=true)

## Installation

1. Use conda virtual environment to install dependencies for lociPARSE. The following command will create a virtual environment named 'lociPARSE'.

```
conda env create -f lociPARSE_environment.yml
```

2. Activate the virtual environment

```
conda activate lociPARSE
```

## Usage

Instructions for running lociPARSE:

1. Put the desired pdb(s) inside the 'Input' folder.

2. Put the PDB ID or list of IDs in the text file named 'input.txt' inside 'Input' folder. See the example in the 'Input' folder.

3. Run
   ```
   chmod a+x lociPARSE.sh && ./lociPARSE.sh
   ```

5. The script will generate features for every ID listed in Input/input.txt and store in individual folder inside 'Feature' folder. Then it will run inference and store predicted molecular-level lDDT (pMoL) and predicted nucleotide-wise lDDT (pNuL) in "score.txt" in individual folder inside 'Prediction' folder.

6. First line in the output "score.txt" shows pMoL score. Each of the subsequent lines sepcify 2 columns: column-1: nucleotide index in PDB and column-2: pNuL score.

## Datasets

- The lists of IDs used in our training set, test sets and validation set used in ablation study are available under [Datasets](https://github.com/Bhattacharya-Lab/lociPARSE/tree/main/Datasets).
- Training set and test set of 30 independent RNAs were taken from [trRosettaRNA](https://yanglab.qd.sdu.edu.cn/trRosettaRNA/benchmark/).
- CASP15 experimental strctures and all submiited predictions were downloaded from [CASP15](https://predictioncenter.org/download_area/CASP15/). 
- 60 non-redundant targets for TS60 validation set were curated from [PDB](https://www.rcsb.org).
