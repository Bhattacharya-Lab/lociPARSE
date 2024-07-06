#!/bin/bash

echo "========Running lociPARSE========="
echo ""
echo "========Generating Features======="

python3 Scripts/Feature.py

echo "Done"

echo "=========Running Inference========"

modelpath=$1

python3 Scripts/prediction.py --modelpath $modelpath

echo ""
echo "=======Inference complete.========"