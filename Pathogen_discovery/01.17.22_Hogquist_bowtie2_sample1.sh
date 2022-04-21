#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu


cd Hogquist_2021/
source activate Hogquist_2021_bowtie2
bowtie2 -p 12 -x GRCm39 -1 1168_1_Qiagen_S1_R1.trim.fastq.gz -2 1168_1_Qiagen_S1_R2.trim.fastq.gz --un-conc-gz 1168_1_Qiagen_S1_host_removed
mv 1168_1_Qiagen_S1_host_removed.1 1168_1_Qiagen_S1_host_removed_R1.fastq.gz
mv 1168_1_Qiagen_S1_host_removed.2 1168_1_Qiagen_S1_host_removed_R2.fastq.gz
