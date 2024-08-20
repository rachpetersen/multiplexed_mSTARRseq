#!/bin/sh

#  mSTARR_trim_map_counts.sh
#  
#
#  Created by Rachel Petersen on 8/15/24.
#  

mSTARR data analysis


##############################trim adapters with cutadapt##############################
cutadapt.sh

#!/bin/bash
#
#SBATCH --mem=20GB
#SBATCH --array=1-494

module load GCC/6.4.0-2.28
module load cutadapt/1.16-Python-3.6.3

fastq_path=/data/lea_lab/petersrm/mSTARR_data/raw
trim_path=/nobackup/lea_lab/petersrm/mSTARR/trimmed

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p Tgen_IDs`

cutadapt -a GATCGGAAGA -A TCTTCCGATC \
-g CTTCCGATCT -G AGATCGGAAG \
-m 20 \
-o ${trim_path}/${sampleID}.1.trimmed.fastq -p ${trim_path}/${sampleID}.2.trimmed.fastq \
${fastq_path}/*${sampleID}*R1*fastq.gz ${fastq_path}/*${sampleID}*R2*fastq.gz

############################## map reads ##############################
map.sh

#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --array=1-494
#SBATCH --mem=32GB
#SBATCH --time=10:00:00

module purge
module load GCC/6.4.0-2.28
module load BWA/0.7.17

fastq_path=/data/lea_lab/petersrm/mSTARR_data/raw
trim_path=/data/lea_lab/petersrm/mSTARR_data/trimmed
genome_path=/data/lea_lab/petersrm/genomes/human
out_path=/nobackup/lea_lab/petersrm/mSTARR/mapped

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p Tgen_IDs`

bwa mem -t 8 ${genome_path}/human_selected.fa \
${trim_path}/${sampleID}.1.trimmed.fastq ${trim_path}/${sampleID}.2.trimmed.fastq \
> ${out_path}/${sampleID}.sam

module purge
module load GCC/10.2.0 SAMtools/1.12

cd mapped

samtools view -b -f 0x2 ${sampleID}.sam | \
samtools sort -n - > ${sampleID}.bam

rm *.sam

## mapped reads
samtools view -F 260 ${sampleID}.bam > ${sampleID}.mapped.bam

## reads uniquely mapped to human
samtools view -h ${sampleID}.mapped.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view > ${sampleID}.human.mapped.uniquely.bam


############################## make genomic windows ##############################

#select only chromosomes (gets rid of unplaced sequences)
cd /data/lea_lab/petersrm/genomes/human
select_chr.sh

#!/bin/bash
module load GCC/10.2.0
module load SAMtools/1.12

samtools faidx hg38.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 \
chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chrX chrY chrMT > human_selected.fa

#Make 400bp windows in genome
makewindows.sh

#!/bin/bash
module load GCC/10.2.0
module load SAMtools/1.12

samtools faidx human_selected.fa
awk -v OFS='\t' {'print $1,$2'} human_selected.fa.fai > human_selected_genomeFile.txt

module purge
module load GCC/8.2.0
module load BEDTools/2.28.0

bedtools makewindows -g human_selected_genomeFile.txt -w 400 > human_selected.400bp.windows.bed

#Make file with 1 unique window per line
awk '{OFS="_";print $1,$2,$3}' human_selected.400bp.windows.bed > human_selected.400bp.windows.txt


############################## make per replicate counts files ##############################
counts_per_rep.sh

#!/bin/bash
#
#SBATCH --array=1-22
#SBATCH --mem=20GB
#SBATCH --time=10:00:00

module purge
module load GCC/10.2.0 SAMtools/1.12
module load Intel/2019.1.144 BEDTools/2.28.0

rep_path=/nobackup/lea_lab/petersrm/mSTARR/replicate_qc

rep=`sed -n ${SLURM_ARRAY_TASK_ID}p replicates`

awk '($3-$2) >= 0 && ($3-$2) <= 2000' $rep_path/${rep}_cat.bed | sed 's/\t*$//' > $rep_path/${rep}_filt.bed

#makes files with 1: chr, 2: start of window, 3: end of window, 4: # of reads in file b that fall in that window, 5: # of bp overlap with that window, 6: # of bases in window (200), 7: percent of bases covered (5/6)
bedtools coverage -a /data/lea_lab/petersrm/genomes/human/human_selected.400bp.windows.bed \
-b $rep_path/${rep}_filt.bed > $rep_path/counts_${rep}_400bpwin.txt

#takes just column 4 with the number of reads falling in that window with no other identifying information
awk '{OFS="\t"; print $4}' $rep_path/counts_${rep}_400bpwin.txt > $rep_path/counts2_${rep}_400bpwin.txt




####################### create counts per rep file for use in R ##############################

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=50GB

module load GCC/10.2.0  OpenMPI/4.0.5
module load R/4.0.5

Rscript make_counts_file.R

## make_counts_file.R

library(data.table)
library(tidyverse)

setwd('/nobackup/lea_lab/petersrm/mSTARR/counts')
b<-read.table('/data/lea_lab/petersrm/genomes/human/human_selected.400bp.windows.txt',header=F,sep='\t')
reps<-read.table('replicates.txt',sep='\t')
reps<-reps$V1
for (p in 1:length(reps)) {
  myfile<-paste0('/nobackup/lea_lab/petersrm/mSTARR/replicate_qc/counts2_',reps[p],'_400bpwin.txt')
  a<-read.table(myfile,header=F,sep='\t')
  colnames(a)<-reps[p]
  b<-cbind(b,a)
}
b2<- b %>%
distinct()
rownames(b2)<-b2[,1]
b2<-b2[,-1]
write.table(b2,file='allcounts_rep_400bpwin.txt',sep='\t',quote=F)


################################ subset DNA and RNA windows to non-zero windows (per replicate) ################################

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=100GB

module load GCC/10.2.0  OpenMPI/4.0.5
module load R/4.0.5

Rscript subset_windows.R

## subset_windows.R

###### subset RNA windows
setwd('/nobackup/lea_lab/petersrm/mSTARR/counts')
info<-read.delim('replicateinfo.txt')

RNAinfo<-subset(info,Sample_type=='RNA')
RNAreps<- RNAinfo$Replicate

Methinfo<- subset(RNAinfo, Meth=="meth")
Methreps<- Methinfo$Replicate

Unmethinfo<- subset(RNAinfo, Meth=="unmeth")
Unmethreps<- Unmethinfo$Replicate

b<-read.table('allcounts_rep_400bpwin.txt, h=T') #filter windows per replicate

RNAcounts<-b[, which(colnames(b) %in% RNAreps)]
meth<-RNAcounts[,which(colnames(RNAcounts) %in% Methreps)]
sham<-RNAcounts[,which(colnames(RNAcounts) %in% Unmethreps)]

meth$zero_counts<-apply(meth,1,function(a) length(which(a==0)))
methfiltered<-subset(meth,zero_counts<3) #Half of reps must be nonzero
methfiltered<-methfiltered[,1:length(methfiltered)-1]

sham$zero_counts<-apply(sham,1,function(a) length(which(a==0)))
shamfiltered<-subset(sham,zero_counts<4) #Half of reps must be nonzero
shamfiltered<- shamfiltered[,1:length(shamfiltered)-1]

dim(shamfiltered)
dim(methfiltered)
write.table(shamfiltered,file='shamrna_halfrepsnonzero_400bpwin.txt',sep='\t')
write.table(methfiltered,file='methrna_halfrepsnonzero_400bpwin.txt',sep='\t')

#If only requiring half nonzero RNA samps in *either* meth or unmeth:
l<-c(rownames(methfiltered),rownames(shamfiltered))
l2<-unique(l)
counts2<-RNAcounts[l2,]
write.table(counts2,file='countsRNArep_halfeithercond_400bpwin.txt',sep='\t')


####### subset DNA windows to non-zero windows
DNAinfo<-subset(info,Sample_type=='DNA')
DNAreps<- DNAinfo$Replicate

Methinfo<- subset(DNAinfo, Meth=="meth")
Methreps<- Methinfo$Replicate

Unmethinfo<- subset(DNAinfo, Meth=="unmeth")
Unmethreps<- Unmethinfo$Replicate

DNAcounts<-b[, which(colnames(b) %in% DNAreps)]
meth<-DNAcounts[,which(colnames(DNAcounts) %in% Methreps)]
sham<-DNAcounts[,which(colnames(DNAcounts) %in% Unmethreps)]

meth$zero_counts<-apply(meth,1,function(a) length(which(a==0)))
methfiltered<-subset(meth,zero_counts<3)
methfiltered<-methfiltered[,1:length(methfiltered)-1]

sham$zero_counts<-apply(sham,1,function(a) length(which(a==0)))
shamfiltered<-subset(sham,zero_counts<4)
shamfiltered<- shamfiltered[,1:length(shamfiltered)-1]

dim(sham)
dim(meth)
write.table(shamfiltered,file='shamdna_halfrepsnonzero_400bpwin.txt',sep='\t')
write.table(methfiltered,file='methdna_halfrepsnonzero_400bpwin.txt',sep='\t')

#If requiring half nonzero DNA samps in *both* meth or unmeth:

counts2<-merge(methfiltered,shamfiltered,by="row.names")
rownames(counts2)<-counts2$Row.names
counts2<-counts2[,-c(1)]
write.table(counts2,file='countsDNArep_halfbothcond_400bpwin.txt',sep='\t')

######## merge DNA and RNA windows

rnaeither<-read.table('countsRNArep_halfeithercond_400bpwin.txt')
#rnaboth<-read.table('countsRNArep_halfbothcond_Msp1.txt')
dna<-read.table('countsDNArep_halfbothcond_400bpwin.txt')

all<-merge(rnaeither,dna,by='row.names')
dim(all)
rownames(all)<-all$Row.names
all<-all[,-c(1)]
write.table(all,file='rnadnaEITHERrep_mergedcounts_400bpwin.txt',sep='\t',quote=F)




## MOVED THESE TO R FOR MODELING ##
