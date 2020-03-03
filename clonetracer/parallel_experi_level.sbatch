#!/bin/bash

#SBATCH -p general
#SBATCH -t 6:00:00
#SBATCH --ntasks=4
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=hhwangbo@unc.edu
#SBATCH --output=parallel_sample_level.log

module add anaconda/2019.10
source /nas/longleaf/apps/anaconda/2019.10/etc/profile.d/conda.sh
source activate clonetracer

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off

#parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog runtask.log --resume"

# this runs the parallel command we want
# in this case, we are running a script named runtask
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}

ls ../../data/clonetracer | grep "DFAM" | parallel -j4 --joblog runtask.log --delay .2 'bash run_experi_level.sh {}'
