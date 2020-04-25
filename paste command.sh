#!/bin/bash
#PBS -S /bin/bash
#PBS -N cut_seq2
#PBS -A kalantry_flux
#PBS -l nodes=1:ppn=4,qos=flux,walltime=08:00:00,mem=16000mb
#PBS -q flux
#PBS -M shriyas@umich.edu
#PBS -m abe
#PBS -V
echo "I ran on:"
cat $PBS_NODEFILE

cd /scratch/kalantry_flux/shriyas/references/strand_filter/
awk 'BEGIN {FS=";"} {OFS="\t"} {print $1}' S14_129_F_intersect > res1
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $11,$20}' S14_129_F_intersect > res2
paste res1 res2 > S14_129_F_cut

awk 'BEGIN {FS=";"} {OFS="\t"} {print $1}' S14_129_B_intersect > res3
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $11,$20}' S14_129_B_intersect > res4
paste res3 res4 > S14_129_B_cut

awk 'BEGIN {FS=";"} {OFS="\t"} {print $1}' S14_jf1_F_intersect > res5
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $11,$20}' S14_jf1_F_intersect > res6
paste res5 res6 > S14_jf1_F_cut

awk 'BEGIN {FS=";"} {OFS="\t"} {print $1}' S14_jf1_B_intersect > res7
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $11,$20}' S14_jf1_B_intersect > res8
paste res7 res8 > S14_jf1_B_cut
