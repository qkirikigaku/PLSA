#!bin/bash
#$ -S /bin/bash

MutData=$1 # MutData: real_sample1 | real_sample2 | sample1

python script/RunMS.py ${MutData} ${SGE_TASK_ID}
