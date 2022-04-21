#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu

cd Hogquist_2021/
source activate Hogquist_2021_bowtie2
bowtie2 -p 12 -x GRCm39 -1 1995_1_Qiagen_S3_R1.trim.fastq.gz -2 1995_1_Qiagen_S3_R2.trim.fastq.gz --un-conc-gz 1995_1_Qiagen_S3_host_removed
mv 1995_1_Qiagen_S3_host_removed.1 1995_1_Qiagen_S3_host_removed_R1.fastq.gz
mv 1995_1_Qiagen_S3_host_removed.2 1995_1_Qiagen_S3_host_removed_R2.fastq.gz
bowtie2 -p 12 -x GRCm39 -1 1995_1_Zymo_S9_R1.trim.fastq.gz -2 1995_1_Zymo_S9_R2.trim.fastq.gz --un-conc-gz 1995_1_Zymo_S9_host_removed
mv 1995_1_Zymo_S9_host_removed.1 1995_1_Zymo_S9_host_removed_R1.fastq.gz
mv 1995_1_Zymo_S9_host_removed.2 1995_1_Zymo_S9_host_removed_R2.fastq.gz

#Sample 5: 1995_1_Qiagen_S3: 1995_1_Qiagen_S3_R1_001.fastq.gz 1995_1_Qiagen_S3_R2_001.fastq.gz \
#Sample 6: 1995_1_Zymo_S9: 1995_1_Zymo_S9_R1_001.fastq.gz 1995_1_Zymo_S9_R2_001.fastq.gz \
