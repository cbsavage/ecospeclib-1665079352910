#!/bin/bash
#SBATCH --account=tkaralidi 
#SBATCH --time=168:00:00
#SBATCH --output=job.%J.out
#SBATCH --job-name=GalileanExoMoons
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


#Select File to run
export file="lava2"
export args=""
#export makefile='makefile'

#Select how logs get stored
mkdir $SLURM_JOB_ID
export debug_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"
export benchmark_logs="$SLURM_JOB_ID/job_$SLURM_JOB_ID.log"

#Load Modules
module load openmpi/openmpi-2.1.1-gcc-7.1.0

# Enter Working Directory
cd $SLURM_SUBMIT_DIR
# Create Log File
echo $SLURM_SUBMIT_DIR
echo "JobID: $SLURM_JOB_ID" >> $debug_logs
echo "Running on $SLURM_JOB_NODELIST" >> $debug_logs
echo "Running on $SLURM_JOB_NNODES nodes." >> $debug_logs
echo "Running on $SLURM_JOB_NPROCS processors." >> $debug_logs
echo  "Current working directory is `pwd`" >> $debug_logs
#$SLURM_NPROCS*$SLURM_NNODES

# Module debugging
module list >> $debug_logs
which mpirun >> $debug_logs

#Start Timestamp
date >> $benchmark_logs
echo "ulimit -l: " >> $benchmark_logs
ulimit -l >> $benchmark_logs
 
# Run job
mpirun ./$file $args
#make
#time -p ./$file
echo "Program is finished with exit code $? at: `date`"
sleep 3

#End Timestamp
echo "Final CPU Hours" >> $benchmark_logs
sacct -u csavage --format=CPUTimeRAW >> $benchmark_logs
date >> $benchmark_logs
echo "ulimit -l" >> $benchmark_logs
ulimit -l >> $benchmark_logs

#Cleanup
#mv job.$SLURM_JOB_ID.err $SLURM_JOB_ID/
#mv job.$SLURM_JOB_ID.out $SLURM_JOB_ID/
#mv fitV2Summary.txt $SLURM_JOB_ID/
