#!/bin/bash

echo "========Predict and Evaluate lociPARSE========="
echo ""


python3 collect.py $1

echo "=========Running Inference========"

modelpath=$2

./lociPARSE_predict.sh $modelpath

cd Evaluation

python3 QA_eval.py $1 1

echo ""
echo "=======Complete.========"