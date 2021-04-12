# ACUgen
Table of Contents
1.	Quick start
2.	Installation
3.	Reference genome
4.	How to use ACUGEN
5.	Example workflows


## Quick Start
•	Install



•	Run ACUGEN

This run should provide the following files:

```
sample.recal.bam
sample.recal.bam.bai
sample_gatk.vcf.gz
sample_freebayes.vcf.gz
sample_bcftools.vcf.gz
sample_merged.btsl.vcf.gz
```

## Installation

ACUGEN is a pipeline that contains more than one genome analysis tools. Before you start using ACUGEN, prerequisites should be provided.
•	Prerequisites

Python 3 (https://www.python.org/)
Docker (https://www.docker.com/)
Burrows-Wheeler Alignment, BWA (http://bio-bwa.sourceforge.net/)
Samtools (http://www.htslib.org/download/)

## Reference genome
We recommend using the GENECODE Human Release 37 (GRCh38.p13) genome for ACUGEN, available below:
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/

GRCh38.p13.genome.fa.gz file must be unziped. Also it must be indexed with BWA and created sequence dictionary file with Samtools before running ACUGEN. Following scripts are down below:

```
bwa index ./GRCh38.p13.genome.fa
samtools dict ./GRCh38.p13.genome.fa
```
## How to use ACUGEN

![image](https://user-images.githubusercontent.com/45030163/114349202-44eda700-9b70-11eb-9e01-ccd61013648c.png)

•	acugen batch – process all FASTQ/BAM files in a directory.
acugen batch fastq – process all paired-end FASTQ files in a directory.
acugen batch bam – process all BAM files in a directory

•	acugen single – process stated FASTQ/BAM files.
acugen single fastq – process stated FASTQ files.
acugen single bam – process stated BAM file.
Only differences of these two function are data formats and how much data you want to run. The both function are provides same outputs. acugen batch finds every file with stated data format and run every file found and, provides FASTQC report, BAM file, individual VCF files of every variant caller in pipeline, merged VCF file for each input. If input file is BAM, no FASTQC report provided. acugen single only runs stated files which user give as input and provide FASTQC report, BAM file, individual VCF files of every variant caller in pipeline, merged VCF file. If input file is BAM, no FASTQC report provided.

 

### acugen batch
acugen batch converts paired-end FASTQ files or BAM files to merged VCF files. Batch function shoud have folder as an input. It reads every input files in the directory and gives different output for each different sample. Paired-end FASTQ files must have “R1” and “R2” in file name after sample name and seperated with “_” from sample name. “bam”, “fastqc”, “individual_vcf” and “merged_vcf” folders will be created to output directory.

acugen batch runs the following flow:

1. Alignment with BWA-MEM
2. Indexing and sorting with Samtools
3. Marking duplicate reads and removing them with Picard
4. Base quality score recalibration with GATK
5. Variant calling with Freebayes, GATK and BCFtools
```
usage:	acugen batch [input format] -i /path/to/input/dir -o /path/to/output/dir

input format	acugen works with FASTQ and BAM files. Input must be “fastq” or “bam”.
-i 		input, path to folder of input files. 
-o 		output, path to folder for output files. If it’s left blank, output files will be created in /acugen/out/
```
### acugen single

acugen single converts paired-end FASTQ files or BAM file to a merged VCF file. Single function requires input file(s) for run. As output “bam”, “fastqc”, “individual_vcf” and “merged_vcf” folders will be created to output directory.
```
Usage: 		acugen single fastq -i <sample_r1.fastq.gz> <sample_r2.fastq.gz> -o /path/to/output/dir
		acugen single bam -i <sample.bam> -o /path/to/output/dir

-i 		input, paired-end FASTQ files or BAM file
-o		output, path to folder for output files. If it’s left blank, output files will be created in /acugen/out/
```


## Example workflows	

Call variants on a single sample paired-end FASTQ files
Use acugen single fastq to create BAM alignment, individual VCFs and merged VCF file from paired-end FASTQ data.

```acugen single fastq -i SAMPLEID_R1.fastq.gz SAMPLEID_R2.fastq.gz -o ./SAMPLEID```

Call variants on a single BAM file
Use acugen single bam to create individual VCFs and merged VCF file from paired-end BAM data.

```acugen single bam -i SAMPLEID.bam -o ./SAMPLEID```

Call variants on different sample paired-end FASTQ files in a directory
Use acugen batch fastq to create BAM alignments, individual VCFs and merged VCF files from paired-end FASTQ files in a directory.

```acugen batch fastq -i ./SAMPLES_FASTQ -o ./SAMPLES_OUT```

Call variants on different sample BAM files in a directory
Use acugen batch bam to create individual VCFs and merged VCF files from BAM files in a directory.

```acugen batch bam -i ./SAMPLES_BAM -o ./SAMPLES_OUT```










