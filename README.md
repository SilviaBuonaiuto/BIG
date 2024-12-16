# The Biorepository and Integrative Genomics (BIG) Initiative

This repository contains the code and workflows used for the analysis and results presented in the paper: 'The Biorepository and Integrative Genomics resource for inclusive genomics: insights from a  diverse pediatric and admixed cohort'. The study investigates the genetic diversity, ancestry stratification, and disease prevalence patterns within the Biorepository and Integrative Genomics (BIG) Initiative, focusing on underrepresented groups in Tennessee, with samples collected specifically from Memphis and its surrounding areas.
The cohort is predominantly pediatric, offering unique insights into early-life health and disease patterns.


# Used tools

A list of tools and software used in this project, along with their versions:

* Python (>= 3.8) : For data processing and analysis script
* R (> 4.0) : For statistical analysis and visualization
* The Burrows-Wheeler Aligner (BWA) MEM : for the alignment of sequence reads to GRCh38 assembly of the human reference genome
* DeepVariant (v0.10.0) : with a custome exome model was used for variant calling
* GLnexus (v1.2.6) : for joint variant calling
* Plink2 : for variants QC and PCA analysis
* bcftools (v1.16) : for variants filtering, QC and vcf manipulation
* Ensembl Variant Effect Predictor (VEP 110): For variants annotation
* SHAPEIT 5 : for Phasing
* ADMIXTURE : for unsupervised cluster analysis
* RFMix (v2.0) : for local and global ancestry inference
* KING : to calculate kinship coefficients
* hap-ibd : to identify IBD segments


