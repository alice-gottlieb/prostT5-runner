#!/bin/bash
# Description: Absolute minimum qsub template.
# You must replace anything in {curly braces} with your own values.
# Lines without any curly braces should be left as is.
# See resource-spec.md for the -l resource values (GPU types, highp, etc.).
#$ -cwd
#$ -o {OUTPUT_LOG_FILE_NAME}
#$ -j y
## h_rt - the maximum runtime of the job. Job is killed after this time
## h_data - the memory *per core*. ie. h_data=4G with one core pulls 4G total.
## e.g. h_rt=23:00:00,h_data=4G
#$ -l {RESOURCES}
#$ -M $USER@ucla.edu
#$ -m ea

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Your Code here
# Always use `uv run python` when running python files
