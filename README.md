
# RNA sequence analysis pipeline using nextflow. 
<img src="https://github.com/dmbala/Nf-RNA-Seq-OARC/blob/master/Fig/dag-flowchart.png" width="450px" height="350px" />

## Files in the repo
 * rna-seq-v5.nf: main pipeline
 * target.csv : design matrix
 * job-submit.s : slurm job description file
 * RNA-seq-nf-timeline.html: Pipeline Execution Summary
 * RNA-seq-nf.html: pipeline MultiQC results

## Major steps 
 * STEP 1 FastQC
 * STEP 2 align with HISAT2
 * STEP 3 RSeQC analysis
 * STEP 4 Mark duplicates
 * STEP 5 Feature counts
 * STEP 6 Merge featurecounts
 * STEP 7 edgeR MDS and heatmap
 * STEP 8 Deseq2 analysis 
 * STEP 9 MultiQC

## Memory bencmarks
The best combinations are
* SBATCH --cpus-per-task=4 and #SBATCH --mem=64GB
or
* SBATCH --cpus-per-task=12 and #SBATCH --mem=192GB

If the alignment or mark duplication step consumes more memory, choose
* SBATCH --cpus-per-task=4 and #SBATCH --mem=192GB


## Single End vs Pairs
params.SingleEnd = "true" for single end
params.SingleEnd = "false" for pairs. 

