#!/bin/bash

#SBATCH -p general
#SBATCH -t 6:00:00
#SBATCH --ntasks=8
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=hhwangbo@unc.edu
#SBATCH --output=parallel_sample_level.log

module add anaconda/2019.10
source /nas/longleaf/apps/anaconda/2019.10/etc/profile.d/conda.sh
source activate clonetracer

# the --exclusive to srun makes srun use distinct CPUs for each job step
# -N1 -n1 allocates a single core to each task
srun="srun --mem=1g -N1 -n1"

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog runtask.log --resume"

# this runs the parallel command we want
# in this case, we are running a script named runtask
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}


#ls ../../data/clonetracer | grep "2020" | $parallel "srun --mem=8g -N1 -n1 bash run_sample_level.sh {}"
#ls ../../data/clonetracer | grep "2020" | $parallel "$srun echo {}"
#for i in $(ls ../../data/clonetracer | grep "2020")
#do
#	#$srun bash run_sample_level.sh $i &i
#	$srun bash -c "cat run_sample_level.sh"
#done

ls ../../data/clonetracer | grep "2020" | parallel -j8 --joblog runtask.log --delay .2 'bash run_sample_level.sh {}'
