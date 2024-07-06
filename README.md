## lociPARSE: a locality-aware invariant point attention model for scoring RNA 3D structures

by Sumit Tarafder and Debswapna Bhattacharya

[[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1)] [[pdf](https://www.biorxiv.org/content/10.1101/2023.11.04.565599v1.full.pdf)]

Codebase for our <ins>loc</ins>ality-aware <ins>i</ins>nvariant <ins>P</ins>oint <ins>A</ins>ttention-based <ins>R</ins>NA <ins>S</ins>cor<ins>E</ins>r (lociPARSE).


<a href="https://zenodo.org/doi/10.5281/zenodo.10369083"><img src="https://zenodo.org/badge/707283184.svg" alt="DOI"></a>


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

Typical installation time on a "normal" desktop computer should take a few minutes in a 64-bit Linux system.

## Usage

Instructions for running lociPARSE:

1. Put the desired pdb(s) inside the 'Input/RNA_pdbs' folder.

2. Put the list of PDB IDs in 'input.txt' inside 'Input' folder. See the example in the 'Input' folder.

3. Run
   ```
   chmod a+x lociPARSE_predict.sh && ./lociPARSE_predict.sh Model/QAmodel_lociPARSE.pt
   ```

5. The script takes the model path as an argument. It will generate features for every ID listed in Input/input.txt and store in individual folder inside 'Feature' folder. Then it will run inference and store predicted molecular-level lDDT (pMoL) and predicted nucleotide-wise lDDT (pNuL) in "score.txt" in individual folder inside 'Prediction' folder.

6. First line in the output "score.txt" shows pMoL score. Each of the subsequent lines sepcify 2 columns: column-1: nucleotide index in PDB and column-2: pNuL score.

Inference time for a typical RNA structure (~70 nucleotides) should take a few seconds.

## Datasets

- The lists of IDs used in our training set, test sets and validation set used in ablation study are available [here](https://zenodo.org/uploads/12669705).
- Training set and test set of 30 independent RNAs were taken from [trRosettaRNA](https://yanglab.qd.sdu.edu.cn/trRosettaRNA/benchmark/).
- CASP15 experimental strctures and all submiited predictions were downloaded from [CASP15](https://predictioncenter.org/download_area/CASP15/). 
- The set of 60 non-redundant RNA targets TS60 for hyperparameter optimization was in-house curated. See (https://doi.org/10.1093/biomethods/bpae047) for more details.

## Training lociPARSE

If you wish to train lociPARSE from scratch on our training set, please follow these steps:

- Download our training dataset train.tar.gz from [here](https://zenodo.org/uploads/12669705) and place it inside Input/Dataset folder.
- Extract the training dataset
  ```
  tar -xzvf train.tar.gz
  ```
-  Run the following command to train our architecture
   ```
   chmod a+x lociPARSE_train.sh && ./lociPARSE_train.sh > log.txt
   ```
   It will take approximately 16 hours to finish feature generation and 50 epochs of training on a single A100 gpu.

- The best model saved on validation loss will be placed inside the Model folder as "QAmodel_retrained.pt". You can use this model to predict as instructed in [Usage](#usage) section.

## Evaluation of lociPARSE

If you want to generate our reported results in the paper from the provided predictions, follow these steps:

-  Extract the provided Evaluation folder which contains all the predictions and ground truths.
   
   ```
   tar -xzvf Evaluation.tar.gz
   ```
-  To generate Tables 1-6, please run the following commands one by one.

   ```
   cd Evaluate
   python3 QA_eval.py Test30_CASP15 0
   python3 QA_eval.py ARES_benchmark2 0
   ```
- You will find the corresponding results inside **Evaluation/Results** folder.
- To generate Supplementary Figures S1-S2, please run the following commands.

  ```
   cd Evaluate
   python3 draw.py
  ``` 
- Generated figures will be inside **Evaluation/Figures** folder.
  
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

-  You will find the corresponding results inside **Evaluation/Results** folder.