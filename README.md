## lociPARSE: a locality-aware invariant point attention model for scoring RNA 3D structures

by Sumit Tarafder and Debswapna Bhattacharya

[[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1)] [[pdf](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1.full.pdf)]

Codebase for our <ins>loc</ins>ality-aware <ins>i</ins>nvariant <ins>P</ins>oint <ins>A</ins>ttention-based <ins>R</ins>NA <ins>S</ins>cor<ins>E</ins>r (lociPARSE).


<a href="https://zenodo.org/doi/10.5281/zenodo.10369083"><img src="https://zenodo.org/badge/707283184.svg" alt="DOI"></a>


![alt text](https://github.com/Bhattacharya-Lab/lociPARSE/blob/main/lociPARSE.png?raw=true)

## Installation
```
pip install lociPARSE
```

Or

```
git clone https://github.com/Bhattacharya-Lab/lociPARSE.git
cd lociPARSE
pip install .
```

Typical installation time should take less than a minute in a 64-bit Linux system.

## Usage

Instructions for running lociPARSE:

```
from lociPARSE import lociparse
lp = lociparse()
score = lp.score("R1108.pdb")
```

Additional functionality

```
score.pMoL.show() #Returns pMoL value
score.pNuL.show() #Returns a list of pNuL values
score.pNuL.show(1) #Returns the pNuL value of 1st nucleotide
score.save("score.txt") #Saves the scores
```

1. Given an RNA pdb "R1108.pdb" as input, lociPARSE predicts both molecular-level lDDT (pMoL) and nucleotide-wise lDDT (pNuL) score.

2. Use show() function to print the pMoL or pNuL values. 

3. Save the output in a provided filename of your choice("score.txt"). First line shows pMoL score. Each of the subsequent lines sepcify 2 columns: column-1: nucleotide index in PDB and column-2: predicted nucleotide-wise lDDT (pNuL) score.

Inference time for a typical RNA structure (~70 nucleotides) should take a few seconds.

## Datasets

- The lists of IDs used in our training set, test sets and validation set used in ablation study are available [here](https://zenodo.org/uploads/12669705).
- Training set and test set of 30 independent RNAs were taken from [trRosettaRNA](https://yanglab.qd.sdu.edu.cn/trRosettaRNA/benchmark/).
- CASP15 experimental strctures and all submiited predictions were downloaded from [CASP15](https://predictioncenter.org/download_area/CASP15/). 
- The set of 60 non-redundant RNA targets TS60 for hyperparameter optimization was in-house curated. See (https://doi.org/10.1093/biomethods/bpae047) for more details.

## Training and evaluation materials

If you want to train or evaluate lociPARSE, please follow these initial steps:

- Download the necessary materials from [here](https://zenodo.org/records/12729167) and place it in the root directory(/lociPARSE)
  ```
  wget https://zenodo.org/records/12729167/files/Materials.tar.gz
  ```

- Extract the Material.tar.gz folder

  ```
  tar -xvzf Materials.tar.gz --strip-components=1
  ```
- Make sure if you have installed appropriate torch version compatible with the CUDA version installed in your machine for GPU training. See here for more (https://pytorch.org/get-started/locally/).

## Training lociPARSE

If you wish to train lociPARSE from scratch on our training set, please follow these steps:

- Download our training dataset Train.tar.gz from [here](https://zenodo.org/uploads/12669705) and place it inside Input/Dataset folder.
- Extract the training dataset
  ```
  tar -xzvf Train.tar.gz
  ```
-  Run the following command to train our architecture
   ```
   chmod a+x lociPARSE_train.sh && ./lociPARSE_train.sh > log.txt
   ```
   It will take approximately 16 hours to finish feature generation and 50 epochs of training on a single A100 gpu.

- The best model saved on validation loss will be placed inside the Model folder as "QAmodel_retrained.pt".

## Evaluation of lociPARSE

If you want to generate our reported results in the paper from the provided predictions, follow these steps:

-  To generate Tables 1-6, please run the following commands one by one.

   ```
   cd Evaluate
   python3 QA_eval.py Test30_CASP15 0
   python3 QA_eval.py ARES_benchmark2 0
   ```
- You will find the corresponding results inside **Evaluate/Results** folder.
- To generate Supplementary Figures S1-S2, please run the following commands.

  ```
   cd Evaluate
   python3 draw.py
  ``` 
- Generated figures will be inside **Evaluate/Figures** folder.
  
If you want to predict the scores by lociPARSE from scratch and re-evaluate, follow these steps: 
 
- Download our test datasets Test.tar.gz and Ares_set.tar.gz from [here](https://zenodo.org/uploads/12669705) and place it inside Input/Dataset folder.

- Extract the folders
  ```
  tar -xzvf Test.tar.gz
  ```
  ```
  tar -xzvf Ares_set.tar.gz
  ```
  
-  To predict and evaluate results on our two test sets Test30 and CASP15 (Tables 1-5), please run the following command.

   ```
   chmod a+x evaluate.sh && ./evaluate.sh Test30_CASP15 Model/QAmodel_lociPARSE.pt
   ``` 

-  To predict and evaluate results on ARES benchmark set-2 (Table 6), please run the following command. [This will be slow due to ~76k models in this test set]

   ```
   chmod a+x evaluate.sh && ./evaluate.sh ARES_benchmark2 Model/QAmodel_Ares_set.pt
   ``` 

-  You will find the corresponding results inside **Evaluate/Results** folder.
