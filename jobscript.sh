#!/bin/bash

#SBATCH -J LPsolver            # Job Name
#SBATCH -A Boltzmann-spectral      # Charges run to <project_account_name>
#SBATCH -o ./output_files/LPsolver.%j.out    # Name of the output and error file (%j expands to jobID) // only symbol % can be recoginized, others say & not
#SBATCH -n 16     # Total number of mpi tasks requested
#SBATCH -N 16       # The number of nodes. n/N tasks are launched on each node. used in conjunction with the -n option (above). Use this option to specify launching less than 16 tasks per node. 
#SBATCH -p development     # Queue (partition) name  -- normal, largemem, development, gpu, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) 
#SBATCH --mail-user=clark.pennie@utexas.edu     #Specify the email address to use for notifications.
#SBATCH --mail-type=ALL       #Specify when user notifications are to be sent.

set -x
export OMP_NUM_THREADS=64
#export I_MPI_DEBUG=2

ibrun tacc_affinity ./bin/LPsolver_2D_Test2D_New2DE_i2toi1.out
# valgrind --leak-check=yes 
#./fpl.out
