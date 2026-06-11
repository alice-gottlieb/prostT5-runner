#!/bin/bash
# Description: Absolutel minimum qsub template.
# You must replace anything in {curly braces} with your own values
# Lines withouth any curly braces should be left as is
#$ -cwd
#$ -o {OUTPUT_LOG_FILE_NAME}
#$ -j y
## h_rt - the maximum runtime of the job. Job is killed after this time
## h_data - the memory *per core*. ie. this script pulls 4G total memory.
#$ -l {RESOURCES} [e.g. h_rt=23:00:00,h_data=4G]
#$ -M $USER@ucla.edu
#$ -m ea

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Your Code here
# Always use `uv run python` when running python files`
