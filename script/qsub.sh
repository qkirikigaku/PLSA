#!bin/bash

qsub -l h_vmem=4G -t 1-50:1 script/RunMS.sh real_sample1
qsub -l h_vmem=8G -t 1-50:1 script/RunMS.sh real_sample2
qsub -l h_vmem=64G -t 1-50:1 script/RunMS.sh sample1

