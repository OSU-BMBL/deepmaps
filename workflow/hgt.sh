#!/usr/bin/bash

cd /home/wan268/hgt/RNA_ATAC
source activate hgt1

# Hyperpapameters
n_batch=50
batch_size=110 #6 10 32
sample_depth=4 #2
n_heads=16
n_layers=2
sample_width=8 #6
rep='T' 
lr=0.2
n_hid=128
epoch=30
reduction='AE'

# File
GAS=/scratch/deepmaps/data/7/GAS.txt
OUT=/scratch/deepmaps/data/7

# Job
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
mkdir $OUT
mkdir $OUT/cell
mkdir $OUT/model
mkdir $OUT/att
mkdir $OUT/gene
python hgt.py --input_dir $GAS --reduction $reduction --n_hid $n_hid --rep $rep --epoch $epoch --lr $lr --n_batch $n_batch --batch_size $batch_size --n_heads $n_heads --n_layers $n_layers --sample_depth $sample_depth --sample_width $sample_width --result_dir $OUT

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
