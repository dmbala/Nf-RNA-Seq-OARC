
/*
===============================================================
 RNA Sequence
===============================================================
 ##  RNA-Seq Analysis Pipeline. Started May 2018.
 ##   https://github.com/dmbala
 ##   Bala Desinghu <dmbala@gmail.com>
 ##   Major steps for RNA-Seq are outlined here at https://github.com/SciLifeLab and more detail about nextflow can be found at https://www.nextflow.io
 ##   ---------------------------------------------------------------
 ##   This pipeline is prepared for 
 ##   the Office of Advanced Research Computing (OARC), Rutgers.  
 ##   The pipeline is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. All required packages are available in OARC cluster. 
 ##   ---------------------------------------------------------------
 #
*/


/* Locations of external programs (rest available via modules) */
Program_Dir = "/projects/oarc/NF-Seq/"

params.reads = "/projects/oarc/NF-Seq/SampleData/chrX_data/samples/*_{1,2}.fastq.gz"
params.SingleEnd=false
params.exp_design_file = "${PWD}/target.csv"

params.genome = "/projects/community/genomics_references/Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
params.gtf = "/projects/community/genomics_references/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf"
params.hisat2_index_base="/projects/oarc/NF-Seq/SampleData/chrX_data/indexes/chrX_tran"
params.outdir = "${PWD}/results"
params.task_cpus = 1
params.do_trimgalore = true

params.do_bbduk = false
params.do_fastp = false
params.do_de_edgeR = false
params.do_de_deseq = true 

params.aligner = "hisat2"
// Strandedness options
forward_stranded = false 
reverse_stranded = false
unstranded = true 


Channel
    .fromFilePairs( params.reads, size: params.SingleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB" }
    .into {read_files_fastqc; read_files_trim_galore; read_files_trim_bbduk; read_files_trim_fastp; ch_print }
ch_print.subscribe {println "read_pairs: $it"}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    $Program_Dir/FastQC/fastqc -q $reads
    """
}

/*
 * STEP 2A - Adaptor trimming and filtering (Trimgalore)
 */


if(params.do_trimgalore){
    process trim_galore {
        tag "Trimgalore"
        module "python/2.7.12"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_galore

        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        file "*trimming_report.txt" into trim_galore_results
        file "*_fastqc.{zip,html}" into trim_galore_fastqc_reports

        script:
        if (params.SingleEnd) {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=$reads out=${name}.clumped.fastq.gz markduplicates subs=0
            $Program_Dir/TrimGalore/trim_galore \\
                --fastqc -q 20 --gzip ${name}.clumped.fastq.gz
            rename _val_1.fq.gz _R1_trimmed.fq.gz *
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/TrimGalore/trim_galore \\
                --paired -q 20 --fastqc \\
                --gzip ${name}.clumped_R1.fastq.gz ${name}.clumped_R2.fastq.gz
            rename _val_1.fq.gz _R1_trimmed.fq.gz *
            rename _val_2.fq.gz _R2_trimmed.fq.gz *
            """
        }

    }
}

/*
 * STEP 2B - Adaptor trimming and filtering (bbduk)
 */

if(params.do_bbduk){
    process trim_bbduk {
        tag "BBDUK"
        module "java/1.8.0_161"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        publishDir "${params.outdir}/trim_bbduk", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".out") > 0) "$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_bbduk

        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        set val(name), file ("*.out") into trim_bbduk_results

        script:
        if (params.SingleEnd) {
            """
            echo "SingleEnd"
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=${reads} out=${name}.clumped.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/bbmap_38.20/bbduk.sh -Xmx8g qin=auto \\
                in=${name}.clumped.fastq.gz \\
                out=${name}_trimmed.fq.gz \\
                ref=$Program_Dir/bbmap_38.20/resources/adapters.fa \\
                stats=${name}.stats.out aqhist=${name}_aqhist.out \\
                ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 ziplevel=9
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/bbmap_38.20/bbduk.sh -Xmx8g qin=auto \\
                in1=${name}.clumped_R1.fastq.gz in2=${name}.clumped_R2.fastq.gz \\
                out1=${name}_R1_trimmed.fq.gz out2=${name}_R2_trimmed.fq.gz \\
                ref=$Program_Dir/bbmap_38.20/resources/adapters.fa \\
                stats=${name}_stats.out aqhist=${name}_aqhist.out \\
                ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 ziplevel=9
            """
        }
    }
check_trim_reads.subscribe {println "trimed read: $it"}
}

/*
 * STEP 2C - Adaptor trimming and filtering (fastp)
 */

if(params.do_fastp){
    process trim_fastp {
        tag "fastp-trimm"
        afterScript = 'rm -rf *.clumped*.fastq.gz'
        module "gcc/5.4"
        publishDir "${params.outdir}/trim_fastp", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_out.html") > 0) "$filename"
                else if (filename.indexOf("_out.json") > 0) "$filename"
                else null
            }

        input:
        set val(name), file(reads) from read_files_trim_fastp


        output:
        set val(name), file ("*trimmed.fq.gz") into trimmed_reads, check_trim_reads, trimmed_reads_fastqc
        set val(name), file ("*_out.*") into trim_fastp_results

        script:
        afterScript = 'rm -rf ${name}.clumped*.fastq.gz'
        if (params.SingleEnd) {
            """
            echo "SingleEnd"
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in=${reads} out=${name}.clumped.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/fastp/fastp -i ${name}.clumped.fastq.gz \\
                -o ${name}_trimmed.fq.gz \\
                -h ${name}_fastp_out.html \\
                -j ${name}_fastp_out.json -z 9
            """
        } else {
            """
            $Program_Dir/bbmap_38.20/clumpify.sh \\
                in1=${reads[0]} in2=${reads[1]}   \\
                out1=${name}.clumped_R1.fastq.gz out2=${name}.clumped_R2.fastq.gz \\
                markduplicates subs=0
            $Program_Dir/fastp/fastp -i ${name}.clumped_R1.fastq.gz \\
                -I ${name}.clumped_R2.fastq.gz -c \\
                -o ${name}_R1_trimmed.fq.gz -O ${name}_R2_trimmed.fq.gz \\
                -h ${name}_fastp_out.html \\
                -j ${name}_fastp_out.json -z 9
            """
        }
    }
check_trim_reads.subscribe {println "trimed read: $it"}
}

/*
 * STEP 3 - Quality control with FastQC after Trim
 */

process trimmed_fastqc {
    tag "trimmed_FastQC"
    publishDir "${params.outdir}/trimmed_fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into trimmed_fastqc_results

    script:
    """
    $Program_Dir/FastQC/fastqc -q $reads
    """
}

if (!params.do_trimgalore && !params.do_bbduk && !params.do_fastp) {
    read_files_trim_galore.into{trimmed_reads; check_trim_reads}
    check_trim_reads.subscribe {println "trim avoided read: $it"}
}


/*
 * STEP 4 - align with HISAT2
 */

hisat2_index_base = params.hisat2_index_base
if(params.aligner == "hisat2"){
    process hisat2Align {
        tag "Hisat2 alignment"
        module "samtools/1.3.1:HISAT2/2.1.0"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: { filename ->
                if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                else null
            }

        input:
        set val(name), file(reads) from trimmed_reads

        output:
        file "${prefix}.bam" into hisat2_bam
        file "${prefix}.hisat2_summary.txt" into alignment_logs

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        def rnastrandness = ''
        if (forward_stranded && !unstranded){
            rnastrandness = params.SingleEnd ? '--rna-strandness F' : '--rna-strandness FR'
        } else if (reverse_stranded && !unstranded){
            rnastrandness = params.SingleEnd ? '--rna-strandness R' : '--rna-strandness RF'
        }
        if (params.SingleEnd) {
            """
            hisat2 -x $hisat2_index_base \\
                   -U $reads \\
                   $rnastrandness \\
                   -p ${params.task_cpus} \\
                   --met-stderr \\
                   --new-summary \\
                   --summary-file ${prefix}.hisat2_summary.txt \\
                   | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
            """
        } else {
            """
            hisat2 -x $hisat2_index_base \\
                   -1 ${reads[0]} \\
                   -2 ${reads[1]} \\
                   $rnastrandness \\
                   --no-mixed \\
                   --no-discordant \\
                   -p ${params.task_cpus} \\
                   --met-stderr \\
                   --new-summary \\
                   --summary-file ${prefix}.hisat2_summary.txt \\
                   | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
            """
        }
    }

    process hisat2_sortOutput {
        tag "HISAT2 sort"
        module "samtools/1.3.1:HISAT2/2.1.0"
        publishDir "${params.outdir}/HISAT2", mode: 'symlink',
            saveAs: { filename -> "aligned_sorted/$filename" }

        input:
        file (reads) from hisat2_bam

        output:
        file "${prefix}.sorted.bam" into bam_count, bam_rseqc, bam_geneBodyCoverage, bam_preseq, bam_markduplicates, bam_ch1

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.bam)?(\.bai)?(\.gz)?$/
        """
        samtools sort \\
            ${reads} \\
            -@ ${params.task_cpus} \\
            -o ${prefix}.sorted.bam
        """
    }
}
bam_ch1.subscribe {println "bam outputs copied: $it"}

/*
 * STEP 5 Build BED12 file
 */


Channel.fromPath(params.gtf)
       .ifEmpty { exit 1, "gtf file not found: ${params.gtf}" }
       .into {gtf_file; gtf_to_bed}

process makeBED12 {
        tag "gen-bed"
        publishDir "${params.outdir}/BedFile", mode: 'copy'

        input:
        file gtf from gtf_to_bed

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: 
        """
        ${Program_Dir}/gtftobed-tool/gtftobed $gtf > ${gtf.baseName}.bed
        """
}


/*
 * STEP 6 - RSeQC analysis
 */

process rseqc {
    tag "${bam_rseqc.baseName - '.sorted'}"
    module "samtools/1.3.1:intel/17.0.2:python/2.7.12"
    module "intel/17.0.4:R-Project/3.4.1"
    
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }

    input:
    file bam_rseqc
    file bed12 from bed_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    def strandRule = ''
    if (forward_stranded && !unstranded){
        strandRule = params.SingleEnd ? '-d ++,--' : '-d 1++,1--,2+-,2-+'
    } else if (reverse_stranded && !unstranded){
        strandRule = params.SingleEnd ? '-d +-,-+' : '-d 1+-,1-+,2++,2--'
    }
    """
    samtools index $bam_rseqc
    set +u
    source ${Program_Dir}/PyProgVirt/rseqc/bin/activate
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    deactivate
    set -u
    """
}


/*
 * Step 6a Rseqc genebody_coverage with subsampling
*/

process genebody_coverage {
    tag "${bam_geneBodyCoverage.baseName - '.sorted'}"
    module "samtools/1.3.1:intel/17.0.2:python/2.7.12"
    module "intel/17.0.4:R-Project/3.4.1"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    input:
    file bam_geneBodyCoverage
    file bed12 from bed_genebody_coverage

    output:
    file "*.{txt,pdf,r}" into genebody_coverage_results

    script:
    """
    cat <(samtools view -H ${bam_geneBodyCoverage}) <(samtools view ${bam_geneBodyCoverage} \\
        | shuf -n 1000000) \\
        | samtools sort - -o ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam
    samtools index ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam
    set +u
    source ${Program_Dir}/PyProgVirt/rseqc/bin/activate
    geneBody_coverage.py \\
        -i ${bam_geneBodyCoverage.baseName}_subsamp_sorted.bam \\
        -o ${bam_geneBodyCoverage.baseName}.rseqc \\
        -r $bed12
    mv log.txt ${bam_geneBodyCoverage.baseName}.rseqc.log.txt
    deactivate
    set -u
    """
}

/*
 * STEP 7 Mark duplicates
 */
process markDuplicates {
    tag "${bam.baseName - '.sorted'}"
    module "java:samtools"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename -> 
            if (filename.indexOf("_metrics.txt") > 0) "$filename"
            else null 
            }


    input:
    file bam from bam_markduplicates

    output:
    file "${bam.baseName}.markDups.bam" into bam_md, bam_featurecounts
    file "${bam.baseName}.markDups_metrics.txt" into picard_results
    file "${bam.baseName}.markDups.bam.bai"

    script:
    """
    java -Xms12G -Xmx12G -XX:ParallelGCThreads=1 -jar ${Program_Dir}/picard-src/bin/picard.jar MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        MAX_RECORDS_IN_RAM=1000000 \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=true \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${bam.baseName}.markDups.bam
    """
}

/*
 * STEP 7 Feature counts
 */

gtf_file = file(params.gtf)
process featureCounts {
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else if (filename.indexOf("_biotype_counts_mqc.txt") > 0) "biotype_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${bam_featurecounts.baseName}_biotype_counts_mqc.txt" into featureCounts_biotype

    script:
    def featureCounts_direction = 0
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    """
    ${Program_Dir}/subread/bin/featureCounts -a $gtf_file -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    ${Program_Dir}/subread/bin/featureCounts -a $gtf_file -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n +3 > ${bam_featurecounts.baseName}_biotype_counts_mqc.txt
    """
}


/*
 * STEP 8 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    module "intel/17.0.2:python/2.7.12"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()

    output:
    file 'merged_gene_counts.txt' into merged_count_edr, merged_count_deseq2

    script:
    """
    ${Program_Dir}/subread/bin/merge_feature_counts.py -o merged_gene_counts.txt -i $input_files
    """
}


/*
 * STEP 9 - edgeR MDS and heatmap
 */

if (params.do_de_edgeR) {
Channel
    .fromPath(params.exp_design_file)
    .ifEmpty { exit 1, "Cannot find any reads matching:" }
    .into {read_file_target_edgeR}

    process edgeR_explore {
        tag "EdgeR exploration analysis on count file"
        module "intel/17.0.4:R-Project/3.4.1"
        publishDir "${params.outdir}/edgeR_exploration", mode: 'copy'

        input:
        file merged_file from merged_count_edr
        file target_input from read_file_target_edgeR

        output:
        file "*.{txt,pdf,csv}" into edgeR_files

        script: // This script is bundled with the pipeline, in NGI-RNAseq/bin/
        """
        Rscript ${Program_Dir}/DE/edgeR_rna-seq_explore.r $merged_file $target_input
        """
    }
}

/*
 * STEP 10 - Deseq2 analysis 
 */

if(params.do_de_deseq) {
Channel
    .fromPath(params.exp_design_file)
    .ifEmpty { exit 1, "Cannot find any reads matching:" }
    .into {read_file_target_deseq2}

    process Deseq2_sample_correlation {
        tag "Deseq2 analysis on count file"
        module "intel/17.0.4:R-Project/3.4.1"
        publishDir "${params.outdir}/deseq2_analysis", mode: 'copy'

        input:
        file merged_file from merged_count_deseq2
        file target_input from read_file_target_deseq2

        output:
        file "*.{txt,pdf,csv}" into deseq2_files

        script: 
        """
        Rscript ${Program_Dir}/DE/deseq2_rna-seq.r $merged_file $target_input
        """
    }
}

/*
 * STEP 11 MultiQC
 */
process multiqc {
    tag "multiqc-results"
    module "intel/17.0.2:python/2.7.12"
    publishDir "${params.outdir}/FullReport-MultiQC", mode: 'copy'

    input:
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    if (params.do_trimgalore) "file ('trim_galore/*') from trim_galore_results.collect()"
    if (params.do_bbduk) "file ('trim_bbduk/*') from trim_bbduk_results.collect()"
    if (params.do_fastp) "file ('trim_fastp/*') from trim_fastp_results.collect()"
    file ('trimmed_fastqc/*') from trimmed_fastqc_results.collect()
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect()
    file ('rseqc/*') from genebody_coverage_results.collect()
    file ('featureCounts/*') from featureCounts_logs.collect()
//toList() as process optional
    if (params.do_de_edgeR) "file ('edgeR_exploration/*') from edgeR_files.toList()" 
    if (params.do_de_deseq) "file ('deseq2_analysis/*') from deseq2_files.toList()"  

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_data

    script:
    rtitle = "RNA-Sequence Analysis on OARC HPC Cluster"
    """
    set +u
    source ${Program_Dir}/PyProgVirt/mqc/bin/activate
    ${Program_Dir}/PyProgVirt/mqc/bin/multiqc -h
    ${Program_Dir}/PyProgVirt/mqc/bin/multiqc -d ${params.outdir} 
    deactivate
    set -u
    """
}


workflow.onComplete {
    println "WORKFLOW SUMMARY"
    println "Pipeline completed at: $workflow.complete"
    println "Duration: $workflow.duration"
    println "WorkDir: $workflow.workDir"
    println "Exit Status: $workflow.exitStatus"
    println ( workflow.complete ? "Okay" : "Oops ..NOT OKAY" )
}


