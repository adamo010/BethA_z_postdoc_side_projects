#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu

cd Hogquist_2021/
source activate Hogquist_2021_bowtie2
bowtie2 -p 12 -x GRCm39 -1 1995_2_Qiagen_S4_R1.trim.fastq.gz -2 1995_2_Qiagen_S4_R2.trim.fastq.gz --un-conc-gz 1995_2_Qiagen_S4_host_removed
mv 1995_2_Qiagen_S4_host_removed.1 1995_2_Qiagen_S4_host_removed_R1.fastq.gz
mv 1995_2_Qiagen_S4_host_removed.2 1995_2_Qiagen_S4_host_removed_R2.fastq.gz
bowtie2 -p 12 -x GRCm39 -1 1995_2_Zymo_S10_R1.trim.fastq.gz -2 1995_2_Zymo_S10_R2.trim.fastq.gz --un-conc-gz 1995_2_Zymo_S10_host_removed
mv 1995_2_Zymo_S10_host_removed.1 1995_2_Zymo_S10_host_removed_R1.fastq.gz
mv 1995_2_Zymo_S10_host_removed.2 1995_2_Zymo_S10_host_removed_R2.fastq.gz

#Sample 7: 1995_2_Qiagen_S4: 1995_2_Qiagen_S4_R1_001.fastq.gz 1995_2_Qiagen_S4_R2_001.fastq.gz \
#Sample 8: 1995_2_Zymo_S10: 1995_2_Zymo_S10_R1_001.fastq.gz 1995_2_Zymo_S10_R2_001.fastq.gz \
