# Variant calling , QC, annotation and phasing

### Alignment
Whole-exome sequence reads where aligned to GRCh38 assembly of the human reference genome
```
bwa mem -t 4 \
GRCh38_full_analysis_set_plus_decoy_hla.fa \
sample.R1.fastq.gz \
sample.R2.fastq.gz \
-R @RG\tID:sample\tPL:ILLUMINA\tPU:puID\tSM:sample\tLB:lbID \
-Y \
-K 100000000  
```
### Variant calling
Single sample variant calling was performed using DeepVariant and joint variant calling using GLnexus workflow

### Variant QC 
QC consists in :
* Excluding variants with more than 5% missing data across all individuals
* Excluding all individuals with more than 10% missing genotypes
* Excluding all variants significantly deviating from Hardy-Weinberg equilibrium (variants with a p-value below 0.0001 for the HWE test will be excluded)
``` 
plink2 \
--vcf BIG.chr${chromosome}.vcf.gz \
--geno 0.05 \
--hwe 0.0001 \
--mind 0.1 \
--recode vcf \
--out BIG.chr${chromosome}.filtered
```
### Variants annotation

### Phasing
* First phase common variants
```
SHAPEIT5_phase_common \
--region chr${chromosome} \
--input BIG.chr${chromosome}.filtered.vcf.gz \
--map 1000Genome.chr${chromosome}.b38.gmap.gz \
--output BIG.chr${chromosome}.phasedCommon.bcf
```
* Phase rare variants using phased common variants as a scaffold
```
while read LINE; do
        CHK=$(echo $LINE | awk '{ print $1; }')
        SRG=$(echo $LINE | awk '{ print $3; }')
        IRG=$(echo $LINE | awk '{ print $4; }')
      	SHAPEIT5_phase_rare \
      	--input BIG.chr${chromosome}.filtered.vcf.gz \
      	--scaffold  BIG.chr${chromosome}.phasedCommon.bcf
      	--map 1000Genome.chr${chromosome}.b38.gmap.gz \
      	--input-region $IRG \ 
      	--scaffold-region $SRG \
      	--output BIG.chr${chromosome}.phased.chunk$CHK.bcf \
      	; done < chr${chromosome}.smallerChunks.txt
```