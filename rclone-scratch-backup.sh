#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/backup-logs/$JOB_ID.out
#$ -j y
#$ -l h_data=4G,h_rt=23:00:00,highp
#$ -pe shared 16
#$ -M $USER@ucla.edu
#$ -m bea
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

/u/home/a/aliceg/bin/rclone copy -P --transfers 16 --checkers 4 --log-file /u/scratch/a/aliceg/logs/backup-logs/scratch-$(date +"%Y-%m-%d_%H-%M-%S").log --log-level INFO /u/scratch/a/aliceg/ ucla-box:research/pellegrini/3di-backups/scratch-$(date +"%Y-%m-%d_%H-%M-%S")