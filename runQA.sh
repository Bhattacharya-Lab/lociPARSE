#!/bin/bash

echo "========Running lociPARSE========="
echo ""
echo "========Generating Features======="

python3 Scripts/Feature.py

echo "Done"

echo "=========Running Inference========"

python3 Scripts/prediction.py
echo ""
echo "=======Inference complete.========"