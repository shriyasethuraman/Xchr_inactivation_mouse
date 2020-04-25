#!/bin/sh
#PBS -S /bin/sh
#PBS -N HTcount
#PBS -A kalantry_flux
#PBS -l nodes=1:ppn=4,qos=flux,walltime=016:00:00,mem=8000mb
#PBS -q flux
#PBS -M shriyas@umich.edu
#PBS -m abe
#PBS -V

echo "I ran on:"
cat $PBS_NODEFILE

module load python
module load med pysam

cd /scratch/kalantry_flux/shriyas/Bioinf545/

mkdir HTcounts
cd HTcounts

python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S1_star/S1_map.bam genes.gtf > HTcounts/count_S1
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S2_star/S2_map.bam genes.gtf > HTcounts/count_S2
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S3_star/S3_map.bam genes.gtf > HTcounts/count_S3
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S4_star/S4_map.bam genes.gtf > HTcounts/count_S4
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S5_star/S5_map.bam genes.gtf > HTcounts/count_S5
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S6_star/S6_map.bam genes.gtf > HTcounts/count_S6
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S7_star/S7_map.bam genes.gtf > HTcounts/count_S7
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S8_star/S8_map.bam genes.gtf > HTcounts/count_S8
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S9_star/S9_map.bam genes.gtf > HTcounts/count_S9
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S10_star/S10_map.bam genes.gtf > HTcounts/count_S10
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S11_star/S11_map.bam genes.gtf > HTcounts/count_S11
python -m HTSeq.scripts.count -f bam -s reverse /scratch/kalantry_flux/shriyas/Bioinf545/star_results/S12_star/S12_map.bam genes.gtf > HTcounts/count_S12


#!/bin/sh
#PBS -S /bin/sh
#PBS -N DEseq
#PBS -A kalantry_flux
#PBS -l nodes=1:ppn=4,qos=flux,walltime=024:00:00,mem=16000mb
#PBS -q flux
#PBS -M shriyas@umich.edu
#PBS -m abe
#PBS -V
echo "I ran on:"
cat $PBS_NODEFILE


cd /scratch/kalantry_flux/shriyas/Michael

module load R


R CMD BATCH --vanilla diff_expr diff_expr_R_out


a<- read.table("count_Ref1_129",header=F,sep="\t",stringsAsFactors=F)
b<- read.table("count_Ref1_jf1",header=F,sep="\t",stringsAsFactors=F)
c<- read.table("count_Ref2_129",header=F,sep="\t",stringsAsFactors=F)
d<- read.table("count_Ref2_jf1",header=F,sep="\t",stringsAsFactors=F)
e<- read.table("count_Ref3_129",header=F,sep="\t",stringsAsFactors=F)
f<- read.table("count_Ref3_jf1",header=F,sep="\t",stringsAsFactors=F)
g<- read.table("count_S1_129",header=F,sep="\t",stringsAsFactors=F)
h<- read.table("count_S1_jf1",header=F,sep="\t",stringsAsFactors=F)
i<- read.table("count_S2_129",header=F,sep="\t",stringsAsFactors=F)
j<- read.table("count_S2_jf1",header=F,sep="\t",stringsAsFactors=F)
k<- read.table("count_S3_129",header=F,sep="\t",stringsAsFactors=F)
l<- read.table("count_S3_jf1",header=F,sep="\t",stringsAsFactors=F)
m<- read.table("count_S4_129",header=F,sep="\t",stringsAsFactors=F)
n<- read.table("count_S4_jf1",header=F,sep="\t",stringsAsFactors=F)


colnames(a) <- c("gene","WT1_129")
colnames(b) <- c("gene","WT1_jf1")
colnames(c) <- c("gene","WT2_129")
colnames(d) <- c("gene","WT2_jf1")
colnames(e) <- c("gene","WT3_129")
colnames(f) <- c("gene","WT3_jf1")
colnames(g) <- c("gene","WT4_129")
colnames(h) <- c("gene","WT4_jf1")
colnames(i) <- c("gene","Mut1_129")
colnames(j) <- c("gene","Mut1_jf1")
colnames(k) <- c("gene","Mut2_129")
colnames(l) <- c("gene","Mut2_jf1")
colnames(m) <- c("gene","Mut3_129")
colnames(n) <- c("gene","Mut3_jf1")


merge1 <- merge(a,b,by=("gene"))
merge2 <- merge(merge1,c,by=("gene"))
merge3 <- merge(merge2,d,by=("gene"))
merge4 <- merge(merge3,e,by=("gene"))
merge5 <- merge(merge4,f,by=("gene"))
merge6 <- merge(merge5,g,by=("gene"))
merge7 <- merge(merge6,h,by=("gene"))
merge8 <- merge(merge7,i,by=("gene"))
merge9 <- merge(merge8,j,by=("gene"))
merge10 <- merge(merge9,k,by=("gene"))
merge11 <- merge(merge10,l,by=("gene"))
merge12 <- merge(merge11,m,by=("gene"))

