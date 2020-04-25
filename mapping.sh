#!/bin/sh
#PBS -S /bin/sh
#PBS -N star_results
#PBS -A kalantry_flux
#PBS -l nodes=2:ppn=4,qos=flux,walltime=6:00:00,mem=32000mb
#PBS -q flux
#PBS -M shriyas@umich.edu
#PBS -m abe
#PBS -V
cd /scratch/kalantry_flux/shriyas/Michael

mkdir star_results
cd star_results

module load med star/2.3.0e

mkdir S12_star

cd S12_star

mkdir S12_129_star

cd S12_129_star

STAR --genomeDir /scratch/kalantry_flux/shriyas/Michael/129_star_ref --readFilesIn /scra\
tch/kalantry_flux/shriyas/Michael/Sample_42814/** /scratch/kalant\
ry_flux/shriyas/Michael/Sample_42814/*** --runThreadN 8 --outFilte\
rMismatchNmax 0 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 3 3 3

module load med samtools
samtools view -bS Aligned.out.sam | samtools sort - S12_129
samtools index S12_129.bam
