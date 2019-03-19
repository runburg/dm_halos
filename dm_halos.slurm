#!/bin/bash
#SBATCH -J dm_halos_hvals.py # Name for your job
#SBATCH -n 1 # Number of tasks when using MPI. Default is 1
#SBATCH -c 1 # Number of cores requested, Default is 1 (total cores requested = tasks x cores)
#SBATCH -N 1 # Number of nodes to spread cores across - default is 1 - if you are not using MPI this should likely be 1
#SBATCH --mem-per-cpu 6400 #the amount of memory per core to request in MB, Default is 3200 MB so be sure to include this (Max for community.q, sb.q and exclusive.q are 6400)
#SBATCH -t 0-10:00 # Runtime in minutes. Default is 10 minutes. The Maximum runtime currently is 72 hours, 4320 minutes - requests over that time will not run
#SBATCH -p community.q# Partition to submit to the standard compute node partition(community.q) or the large memory node partition(lm.q)
#SBATCH -o hvals.out # Standard out goes to this file
#SBATCH -e hvals.err # Standard err goes to this file
#SBATCH --mail-user runburg@hawaii.edu # this is the email you wish to be notified at
#SBATCH --mail-type ALL # this specifies what events you should get an email about ALL will alert you of job beginning, completion, failure etc
source ~/.bash_profile #if you want to use modules or need environment variables use this if your shell is bash to load those
module load lang/Python/3.5.1/python  #if you want to load a module use this and the above common with the module name you wanted loaded.
python3 h_values.py
