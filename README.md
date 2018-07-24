
# RNA sequence analysis pipeline using nextflow. 
<img src="https://github.com/dmbala/Nf-RNA-Seq-OARC/blob/master/Fig/dag-flowchart.png" width="450px" height="350px" />

## Files in the repo
 * rna-seq-v5.nf: main pipeline
 * target.csv : design matrix
 * job-submit.s : slurm job description file
 * RNA-seq-nf-timeline.html: Pipeline Execution Summary
 * RNA-seq-nf.html: pipeline MultiQC results

## Major steps 
 * STEP 1 Quality control check with FastQC
 * STEP 2 Remove adaptors with trimgalore
 * STEP 3 Align with HISAT2
 * STEP 4 Generate BED from gtf file
 * STEP 5 Perform RSeQC analysis after alighment
 * STEP 6 Mark duplicates
 * STEP 7 Feature counts
 * STEP 8 Merge featurecounts
 * STEP 9 edgeR MDS and heatmap
 * STEP 10 Deseq2 analysis 
 * STEP 11 MultiQC

## Memory bencmarks
The best combinations are
* SBATCH --cpus-per-task=4 and #SBATCH --mem=64GB
or
* SBATCH --cpus-per-task=12 and #SBATCH --mem=192GB

If the alignment or mark duplication step consumes more memory, choose
* SBATCH --cpus-per-task=4 and #SBATCH --mem=192GB


## Single End vs Pairs
For single end reads
```
params.SingleEnd = "true"
```
For pair reads
```
params.SingleEnd = "false" 
```
