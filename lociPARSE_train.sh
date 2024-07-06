#!/bin/bash

STARTTIME="$(date -u +%s)"

echo "========Training lociPARSE network========="
echo ""
echo "========Generating Features======="

K=20 #Number of nearest neighbours for IPA

python3 Scripts/Feature.py --filename Input/Dataset/Train/all.txt --num_neighbours $K --mode "train"

echo "Done"

echo "=========Training starts=========="

#==============Hyper Parameters=============
num_layers=4
hidden_dim=128
gpunum=2
epochs=50
lr_rate=1e-4
wdecay=5e-6
droprate=0.1
SEED=1000


python3 -u Scripts/train.py --epochs $epochs --lr $lr_rate --n_layers $num_layers --nf $hidden_dim --weight_decay $wdecay --droprate $droprate --gpunum $gpunum --numneighbour $K --SEED $SEED

echo ""
echo "=======Training complete.========="

ENDTIME="$(date -u +%s)"

echo "Total time taken to train = $(( ($ENDTIME - $STARTTIME) / 3600 )) hours."
