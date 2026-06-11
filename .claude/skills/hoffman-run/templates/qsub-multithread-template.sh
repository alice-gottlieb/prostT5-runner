#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
## h_rt - the maximum runtime of the job. Job is killed after this time
## h_data - the memory *per core*. ie. this script pulls 2G * 4 = 8G total memory.
#$ -l h_rt=1:00:00,h_data=2G
## pe - parallel environment. 
## shared - allows multiple jobs to share the same node and share memory. 
## Number specified is the number of cores.
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
## Send email at the beginning (b), end (e), and if the job is aborted (a)
## Can remove 'b' if you don't want to receive an email when the job starts, 
## 'e' if you don't want to receive an email when the job ends, 
## and 'a' if you don't want to receive an email if the job is aborted
#$ -m bea

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Your Code here