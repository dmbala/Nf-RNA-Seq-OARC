#!/bin/bash
#SBATCH --job-name=check-seq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=6:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
module use /projects/community/modulefiles
module load nextflow 
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir
#srun nextflow run rna-seq-v5.nf --task_cpus=$SLURM_CPUS_PER_TASK --SingleEnd="true"  -w $NF_Work_Dir -with-trace -with-report RNA-seq-nf.html  -with-timeline RNA-seq-nf-timeline.html -resume 
#srun nextflow run rna-seq-v5.nf --SingleEnd="false"  --reads=/projects/oarc/NF-Seq/SampleData/chrX_data/samples-orig/*.fastq.gz" -w $NF_Work_Dir -with-trace -with-report RNA-seq-nf.html  -with-timeline RNA-seq-nf-timeline.html -resume -with-dag dag-flowchart.html
srun nextflow run rna-seq-v5.nf --SingleEnd="false"  -w $NF_Work_Dir -resume -with-trace 

if  [ $?  -eq  0 ];
then
    echo " Pipeline completed. Removing WorkDir files"
    rm -rf $NF_Work_Dir
else
    echo "Pipeline not completed." 
fi


