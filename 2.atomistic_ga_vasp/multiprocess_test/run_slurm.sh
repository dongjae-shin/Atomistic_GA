#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20     # Cores per node
#SBATCH --partition=g3           # Queue Name g3
##
#SBATCH --job-name="GA_test"
#SBATCH --time=7-12:30          # Runtime limit: Day-HH:MM
#SBATCH -o STDOUT.%N.%j.out         # STDOUT, %N : nodename, %j : JobID
#SBATCH -e STDERR.%N.%j.err         # STDERR, %N : nodename, %j : JobID


## HPC ENVIRONMENT DON'T REMOVE THIS PART
. /etc/profile.d/TMI.sh


#mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS  /TGM/Apps/VASP/bin/5.3.5/NORMAL/vasp.5.3.5_31MAR2014_GRP7_NORMAL_VTST.x
#mpiexec.hydra -genv I_MPI_DEBUG 5 -np $SLURM_NTASKS /home/djshin/research/python_tests/wait/wait_test.py > stdout.txt
#python -u /home/djshin/research/python_tests/wait/wait_test.py > stdout.txt

python -u /home/djshin/research/1.181115_GA/python_tests/multiprocess_test/multiprocessing_test.py > stdout.txt
