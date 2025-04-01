#!/bin/sh

#  GATK_joint_genotyping_ASE.sh
#  
#
#  Created by Rachel Petersen on 4/1/25.
#  


#### GATK snp calling

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --array=1

module purge
module load GCC/8.2.0 SAMtools/1.9 picard/2.18.27

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p human_ids`

in_bam=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/18Jul_IDchanges/${sampleID}_dna_p1.bam
in_genome=/data/lea_lab/petersrm/genomes/human/human_selected.fa
in_variants=/data/lea_lab/shared/1KG_mSTARR_variants/1000G_phase1.snps.high_confidence.hg38.vcf.gz

out_bam_sorted=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/18Jul_IDchanges/$sampleID_dna_p1.sorted.bam
out_bam_duplicates=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/$sampleID.dna.p1.dup.bam
out_txt_dupmetrics=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/$sampleID.dna.p1.dup.txt

out_bam_readgroups=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/$sampleID.dna.p1.reg.bam
out_table_recalibr=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/$sampleID.dna.p1.rbqs.table
out_bam_recalibrat=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/$sampleID.dna.p1.rbqs.bam
out_g_vcf=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/$sampleID.p1.g.vcf

# Step   I: Identify duplicated reads with sorted bam or sam files
# Step  II: Recalibrate base quality scores
# Step III: Calling GATK variants and export VCF

# Step I
echo 'samtools sorting...' # sort bam file

samtools sort -m 3G -o $out_bam_sorted $in_bam

echo 'picard duplications...' # run picard to mark duplicates

java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$out_bam_sorted O=$out_bam_duplicates M=$out_txt_dupmetrics

# Step II
echo 'adding read groups...' ##adds a header

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$out_bam_duplicates O=$out_bam_readgroups SO=coordinate RGLB=$sampleID RGPL=illumina RGPU=mSTARR RGSM=$sampleID

rm -f $out_bam_duplicates
rm -f $out_txt_dupmetrics

echo 'pre-indexing...'

samtools index $out_bam_readgroups

echo 'recalibrating base quality scores...' # sort bam file

gatk BaseRecalibrator -I $out_bam_readgroups -R $in_genome \
--known-sites $in_variants -O $out_table_recalibr

gatk ApplyBQSR -I $out_bam_readgroups -bqsr $out_table_recalibr -O $out_bam_recalibrat

rm -f $out_bam_readgroups

echo 'post-indexing...'

samtools index $out_bam_recalibrat

# Step III
echo 'calling GATK variants as VCF and zip...' # sort bam file

gatk HaplotypeCaller -R $in_genome -I $out_bam_recalibrat -O $out_g_vcf -ERC GVCF --max-alternate-alleles 2

module purge
module load Intel/2016.3.210 tabix/0.2.6

bgzip -f $out_g_vcf
tabix -f -p vcf $out_g_vcf.gz

rm $out_bam_recalibrat


##### joint genotyping

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem=10GB

module purge
module load GCCcore/.8.2.0
module --ignore_cache load GATK

chr=`sed -n ${SLURM_ARRAY_TASK_ID}p chrom_intervals`

chrom=$chr
in_genome=/data/lea_lab/petersrm/genomes/human/human_selected.fa
in_database=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/joint_geno/$chrom
in_samples=/nobackup/lea_lab/petersrm/mSTARR/sample-name-map-mergedlibs.txt
out_vcf_gz=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/joint_geno/$chrom.vcf.gz
tmp_dir=/nobackup/lea_lab/petersrm/mSTARR/tmp_files

# merge GVCFs
gatk --java-options "-Xmx80g -Xms80g" GenomicsDBImport \
            --sample-name-map $in_samples \
            --genomicsdb-workspace-path $in_database \
            --tmp-dir=$tmp_dir \
            -L $chrom \
            --batch-size 40 --reader-threads 16

# call genotypes

echo 'calling joint genotypes...'

echo $chrom > temp.${chrom}.list

gatk --java-options "-Xmx80g" GenotypeGVCFs \
            -R $in_genome \
            -L temp.${chrom}.list \
            -V gendb://$in_database \
            -O $out_vcf_gz \
            --tmp-dir=$tmp_dir --max-alternate-alleles 2


#### concatenate vcf files of each chromosome and intersect joint genotyping vcf with 1000G vcf

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=50GB

module purge
module load GCC/10.2.0 BCFtools/1.16

vcf_path=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/joint_geno/18Jul_IDchanges

bcftools concat -o $vcf_path/allCHR.vcf.gz $vcf_path/chr*.vcf.gz

tabix -f -p vcf $vcf_path/allCHR.vcf.gz

bcftools isec -p isec_output -Oz $vcf_path/allCHR.vcf.gz /data/lea_lab/shared/1KG_mSTARR_variants/1000G_phase1.snps.high_confidence.hg38.vcf.gz

#0000.vcf.gz is variants unique to 1.vcf.gz
#0001.vcf.gz is variants unique to 2.vcf.gz
#0002.vcf.gz is variants shared by 1.vcf.gz and 2.vcf.gz as represented in 1.vcf.gz
#0003.vcf.gz is variants shared by 1.vcf.gz and 2.vcf.gz as represented in 2.vcf.gz


# count number of lines that do not start with #
zcat 0002.vcf.gz | grep "^[^#]" | wc -l
#737,394


######## merge bam files within a replicate

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --array=1

repID=`sed -n ${SLURM_ARRAY_TASK_ID}p replicates`

module load GCC/8.2.0 SAMtools/1.9

samtools merge /nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/${repID}.bam -b /nobackup/lea_lab/petersrm/mSTARR/replicate_qc/${repID}_bams


##### Allele specific expression

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=50GB
#SBATCH --array=1

repID=`sed -n ${SLURM_ARRAY_TASK_ID}p replicates`

in_gatk=/nobackup/lea_lab/petersrm/mSTARR/gatk-4.1.4.0
in_bam=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/${repID}.bam
out_bam_sorted=/nobackup/lea_lab/petersrm/mSTARR/allele_based_analyses/${repID}.sorted.bam
in_genome=/data/lea_lab/petersrm/genomes/human/human_selected.fa
in_vcf=/nobackup/lea_lab/petersrm/mSTARR/isec_output/0002.vcf.gz
out_bam_readgroups=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/${repID}.rg.bam
out_ASE=/nobackup/lea_lab/petersrm/mSTARR/variant_calling/${repID}.ase
windows=/nobackup/lea_lab/petersrm/mSTARR/OneHapWindows.bed

module purge
module load GCC/8.2.0 picard/2.18.27 SAMtools/1.9

samtools sort $in_bam > \
$out_bam_sorted
samtools index $out_bam_sorted

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=$out_bam_sorted \
O=$out_bam_readgroups \
SO=coordinate \
RGLB=rep1 \
RGPL=illumina \
RGPU=mSTARR \
RGSM=$repID

samtools index $out_bam_readgroups

$in_gatk/gatk ASEReadCounter -I $out_bam_readgroups \
-R $in_genome \
-V $in_vcf \
-DF NotDuplicateReadFilter \
-L $windows \
-O $out_ASE

