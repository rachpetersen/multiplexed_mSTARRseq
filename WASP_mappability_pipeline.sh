#!/bin/sh

#  WASP_mappability_pipeline.sh
#  
#
#  Created by Rachel Petersen on 4/1/25.
#  

#Step 1: make SNP text files

INPUT_DIR=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/joint_geno/18Jul_IDchanges
OUTPUT_DIR=/nobackup/lea_lab/petersrm/mSTARR/WASP/output_snp_dir

for FILE in $INPUT_DIR/*vcf.gz; do
     echo $FILE >&2
     CHR=`echo $FILE | sed -n 's/^.*\(chr[0-9A-Z]*\).*.vcf.gz$/\1/p'`
     echo $CHR >&2
     OUTPUT_FILE=$OUTPUT_DIR/$CHR.snps.txt.gz
     gunzip -c $FILE | egrep -v "^#" | awk '{print $2,$4,$5}' | gzip > $OUTPUT_FILE
done

#Step 2: sort and index existing mapped files

#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --array=2-494
#SBATCH --mem=2GB
#SBATCH --time=00:10:00

module load GCC/11.3.0 SAMtools/1.18

in_path=/nobackup/lea_lab/petersrm/mSTARR/mapped/uniquely_mapped
out_path=/nobackup/lea_lab/petersrm/mSTARR/mapped/uniquely_mapped

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p Tgen_IDs`


samtools sort -o ${out_path}/${sampleID}.human.mapped.uniquely.sort.bam ${in_path}/${sampleID}.human.mapped.uniquely.bam
samtools index ${out_path}/${sampleID}.human.mapped.uniquely.sort.bam


#Step 3: find_intersecting_snps.sh

#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --array=2
#SBATCH --mem=10GB
#SBATCH --time=00:20:00

module load Anaconda3/2023.03-1
module load GCC/11.3.0 Pysam/0.19.1

module load GCC/8.2.0 Pysam/0.16.0.1
module load GCC/6.4.0-2.28 Pysam/0.15.2-Python-3.6.3


sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /nobackup/lea_lab/petersrm/mSTARR/Tgen_IDs`


python mapping/find_intersecting_snps.py \
          --is_paired_end \
          --is_sorted \
          --output_dir find_intersecting_snps \
          --snp_dir /nobackup/lea_lab/petersrm/mSTARR/WASP/output_snp_dir \
          /nobackup/lea_lab/petersrm/mSTARR/mapped/uniquely_mapped/L51775.human.mapped.uniquely.sort.bam

#Step 4: Remap using WASP_remap.sh

#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --array=1
#SBATCH --mem=32GB
#SBATCH --time=2:00:00

module purge
module load StdEnv/2023
module load bwa/0.7.18

fastq_path=/nobackup/lea_lab/petersrm/mSTARR/WASP/find_intersecting_snps
genome_path=/data/lea_lab/shared/genomes/hg38
out_path=/nobackup/lea_lab/petersrm/mSTARR/WASP/remap

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p Tgen_IDs`

#bwa index ${genome_path}/hg38_ALLchr.fa

bwa mem -t 8 ${genome_path}/hg38_ALLchr.fa \
${fastq_path}/${sampleID}.human.mapped.uniquely.sort.remap.fq1.gz ${fastq_path}/${sampleID}.human.mapped.uniquely.sort.remap.fq2.gz \
> ${out_path}/${sampleID}.sam

module load gcc/12.3 samtools/1.20

cd /nobackup/lea_lab/petersrm/mSTARR/WASP/remap

samtools view -b -f 0x2 ${sampleID}.sam | \
samtools sort -o ${sampleID}.remap.sorted.bam
samtools index ${sampleID}.remap.sorted.bam


#Step 5: filter remapped reads

python mapping/filter_remapped_reads.py \
       find_intersecting_snps/${sampleID}.human.mapped.uniquely.sort.to.remap.bam \
       remap/${sampleID}.remap.sorted.bam \
       filter_remapped_reads/${sampleID}.keep.bam
       
       
 #Step 6: merge
 
#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --array=1-494
#SBATCH --mem=5GB
#SBATCH --time=00:30:00
#SBATCH --constraint=zen
 
module load StdEnv/2023 gcc/12.3 samtools/1.20
 
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /nobackup/lea_lab/petersrm/mSTARR/Tgen_IDs`
 
samtools merge merge/${sampleID}.keep.merge.bam \
              filter_remapped_reads/${sampleID}.keep.bam  \
              find_intersecting_snps/${sampleID}.human.mapped.uniquely.sort.keep.bam
samtools sort -o  merge/${sampleID}.keep.merge.sort.bam \
              merge/${sampleID}.keep.merge.bam
samtools index merge/${sampleID}.keep.merge.sort.bam

