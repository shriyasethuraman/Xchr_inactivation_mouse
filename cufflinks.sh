#!/bin/sh
#PBS -S /bin/sh
#PBS -N cufflinks_rep1
#PBS -A kalantry_flux
#PBS -l nodes=1:ppn=4,qos=flux,walltime=024:00:00,mem=16000mb
#PBS -q flux
#PBS -M shriyas@umich.edu
#PBS -m abe
#PBS -V
echo "I ran on:"
cat $PBS_NODEFILE

module load boost
module load med samtools/0.1.18
module load med cufflinks/2.2.1

cd /scratch/kalantry_flux/shriyas/references/star_transcriptome

cufflinks -o ./cufflinks_rep1 -p 4 --library-type fr-firststrand --min-frags-per-transfrag 5 /scratch/kalantry_flux/shriyas/references/star_results/rep1_star/rep1_129.bam