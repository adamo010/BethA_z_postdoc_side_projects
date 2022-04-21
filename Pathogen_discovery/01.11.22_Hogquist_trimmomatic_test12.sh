#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu

cd Hogquist_2021/
source activate Hogquist_2021
trimmomatic PE Extraction_BLANK_2_Zymo_S12_R1_001.fastq.gz Extraction_BLANK_2_Zymo_S12_R2_001.fastq.gz \
                Extraction_BLANK_2_Zymo_S12_R1.trim.fastq.gz Extraction_BLANK_2_Zymo_S12_R1.untrim.fastq.gz \
                Extraction_BLANK_2_Zymo_S12_R2.trim.fastq.gz Extraction_BLANK_2_Zymo_S12_R2.untrim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
