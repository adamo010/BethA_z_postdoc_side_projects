#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu

cd Hogquist_2021/
source activate Hogquist_2021_bowtie2
bowtie2 -p 12 -x GRCm39 -1 Extraction_BLANK_1_Qiagen_S5_R1.trim.fastq.gz -2 Extraction_BLANK_1_Qiagen_S5_R2.trim.fastq.gz --un-conc-gz Extraction_BLANK_1_Qiagen_S5_host_removed
mv Extraction_BLANK_1_Qiagen_S5_host_removed.1 Extraction_BLANK_1_Qiagen_S5_host_removed_R1.fastq.gz
mv Extraction_BLANK_1_Qiagen_S5_host_removed.2 Extraction_BLANK_1_Qiagen_S5_host_removed_R2.fastq.gz
bowtie2 -p 12 -x GRCm39 -1 Extraction_BLANK_1_Zymo_S11_R1.trim.fastq.gz -2 Extraction_BLANK_1_Zymo_S11_R2.trim.fastq.gz --un-conc-gz Extraction_BLANK_1_Zymo_S11_host_removed
mv Extraction_BLANK_1_Zymo_S11_host_removed.1 Extraction_BLANK_1_Zymo_S11_host_removed_R1.fastq.gz
mv Extraction_BLANK_1_Zymo_S11_host_removed.2 Extraction_BLANK_1_Zymo_S11_host_removed_R2.fastq.gz

#Sample 9: Extraction_BLANK_1_Qiagen_S5: Extraction_BLANK_1_Qiagen_S5_R1_001.fastq.gz Extraction_BLANK_1_Qiagen_S5_R2_001.fastq.gz \
#Sample 10: Extraction_BLANK_1_Zymo_S11: Extraction_BLANK_1_Zymo_S11_R1_001.fastq.gz Extraction_BLANK_1_Zymo_S11_R2_001.fastq.gz \
