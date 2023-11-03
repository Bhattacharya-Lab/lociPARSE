# lociPARSE: a locality-aware invariant point attention model for scoring RNA 3D structures

by Sumit Tarafder and Debswapna Bhattacharya.

Submitted, 2023.


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

