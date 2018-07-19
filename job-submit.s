#!/bin/bash
#SBATCH --job-name=check-seq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --time=24:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
module use /projects/community/modulefiles
module load nextflow 
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir
#srun nextflow run rna-seq-v5.nf --task_cpus=$SLURM_CPUS_PER_TASK --SingleEnd="true"  -w $NF_Work_Dir -with-trace -with-report RNA-seq-nf.html  -with-timeline RNA-seq-nf-timeline.html -resume 
srun nextflow run rna-seq-v5.nf --SingleEnd="true"  -w $NF_Work_Dir -with-trace -with-report RNA-seq-nf.html  -with-timeline RNA-seq-nf-timeline.html -resume 

