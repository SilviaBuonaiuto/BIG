# LocalAncestry
This repo contains a pipeline for local ancestry estimation using rfmix2 and 1000GP  (plus HGDP, optional) data


# Local ancestry estimation
Firstly, for natural selection analyses, local ancestry estimation should be done for admixed individuals. The following roadmap is based on this [protocol](https://star-protocols.cell.com/protocols/1034), but adapted for refMix2.


### Pipeline:

- **Local Ancestry Estimation**:
    1. **Data Preparation**:
        - Low-quality SNPs or individuals with low call rates are removed using [PLINK](http://zzz.bwh.harvard.edu/plink/).
          
    2. **Reference Populations**:
        - Identify reference populations that represent the sources of ancestry in the study population. Obtain samples from 1000 Genomes Project, The Human Genome Diversity Project, or other. Particularly, joint call is a very useful reference: [Data](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg-tutorials) and [Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9900804/)
        - Reference genome: GRCh38/hg38.
        - Perform a QC on reference populations. Firstly, with unsupervised clustering (ADMIXTURE) select those individuals with more than 99% of a single component. 
          
    3. **Phasing**:
        - Phase the genotypes [SHAPEIT](https://odelaneau.github.io/shapeit4/).
        - Software: High-performance computing (HPC) environment.
          
    4. **Local Ancestry Estimation**:
        - [RFMix](https://github.com/slowkoni/rfmix) to estimate local ancestry in the population. It should be done in HPC environment.
          
In this case the reference selected population are: AFR reference populations YRI MSL GWD;  EUR ref population CEU TSI; AMR refrence population PEL ;   EAS reference CHB [www.internationalgenome.org]


### QC:

```bash
plink --geno 0.05 # (filtered out all variants with missing call rates exceeding the provided value to be removed)
plink --hwe 0.0001 #(filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold)
plink --mind 0.1 #(filters out all samples with missing call rates exceeding the provided value to be removed)
```

### **SHAPEIT5**:
```
Phase common variants (MAF >= 0.1%)
Phase rare variants (MAF <= 0.1%) using phased common variants as scaffold
```

### **Obtaining reference populations data**:

For this point, we will download phased vcf from 1000GP from the following link: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/. 
The following code could be executed for downloading the data:

For this pilot analysis, we will focus on chr15:

```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

# Once downloaded, we can name both files as AllPop.vcf.gz and AllPop.vcf.gz.tbi. Then, we can take a set of samples as reference

mv 1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz AllPop.vcf.gz
mv 1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi AllPop.vcf.gz.tbi
```

# Pilot (for testing the software rapidly and in a more controlled sample)
Estimated time, 64gb ram, 20 cores: 35 mins


```bash
bcftools view -s NA19836,NA19901,NA19913,NA19918,NA19920,NA19982,NA20129,NA20282,NA19704,NA19711,NA11829,NA12342,NA11831,NA12347,NA11843,NA12400,NA11881,NA11893,NA12739,NA11918 AllPop.vcf.gz -Oz -o Reference.vcf.gz
```

Then we index the data:
```bash
tabix -p vcf Reference.vcf.gz
```

NOTE: this samples are also listed in the sample_pilot.map file. All files should be in a phased .vcf.gz format.

Just as an example, if we want to test the pipeline with online data: Lets supose that our Query is a person from Puerto Rico:

```bash
bcftools view -s HG00638 AllPop.vcf.gz -Oz -o Query.vcf.gz
tabix -p vcf Query.vcf.gz
```


If everything works well, we are able to start running RFMix2. 


## **Input preparation**:
We need some inputs for running RFmix2.

1. **Query.vcf.gz**: 
   - This is the cohort, that could be phased independently.

2. **Reference.vcf.gz**: 
   - This is the phased vcf with data from the two reference populations, here ASW and CEU. previously downloaded.
     
3. **sample map**: 
   - This file must contain two columns \tab delimited. First column indicates the sample name (it should be not repeted within and between reference and query populations). It is already available for the pilot study. It is named sample_pilot.map.

4. **Genetic map**: Also we need the genetic map, that could be obtained from the following link [genetic map hg18](https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg18.txt.gz). But, if it is used directly, will give an strange error: "Stop: chr is not striclty increasing". This could be repeared extracting columns 1, 3 and 4 from the genetic map and changing the number of the chromosome, for example 15, by chr15. It is already available on this link: https://drive.google.com/file/d/1_1dVIi2de8ZaH4P0j-eFj2jxAMqJAiJN/view?usp=sharing

## **Install RFmix**:
Just download the repository from https://github.com/slowkoni/rfmix
```bash
git clone https://github.com/slowkoni/rfmix
```

and then run:

```bash
cd rfmix \
autoreconf --force --install # creates the configure script and all its dependencies \
./configure                  # generates the Makefile \
make \
cd ..
```

It will produce rfmix executable.

## **Infering local ancestry with RFmix**:

Here we will put the code for a pilot study, just focusing on chrm15 (it is short), and a low number of unrelated individuals in the reference (5 CEU and 5 ASW). The command for running RFmix2 is the following:

```bash
./rfmix/rfmix -f Query.vcf.gz -r Reference.vcf.gz -m sample_pilot.map -g genetic_map_hg18_rfmix_modified.txt -o local15_PR --chromosome=chr15 --n-threads=20
```

It should run relatively fast (around 20 minutes in a 64gb ram, 20 cores used). The following files will be obtained:

1. local15_PR.msp.tsv: Important for karyogram
2. local15_PR.rfmix.Q: ancestry proportion
3. local15_PR.sis.tsv and local15_PR.fb.tsv: other information related to specific snps ancestry.

# Running a full analysis (1 chromosome - two source populations)
Estimated time, 64GB RAM, 20 cores: 45 mins per chromosome 15.

Slight modifications have been made to run the analysis with all the samples. In these cases, we continue using Puerto Rico as the target population.

NOTE: We continue with chromosome 15. We can easily run all together (but it will probably be useful to prioritize, as there is too much data for analysis).

We download (or use, if previously downloaded) the same data obtained with wget (first command lines of this README file). bcftools view -S command allows extracting a set of specified samples (listed in reference_list.txt) from a .vcf.gz file (AllPop.vcf.gz, downloaded with wget previously, see Obtaining reference populations data section).

```bash
bcftools view -S query_list.txt AllPop.vcf.gz -Oz -o Query_full.vcf.gz --force-samples
bcftools view -S reference_list.txt AllPop.vcf.gz -Oz -o Reference_full.vcf.gz --force-samples
tabix -p vcf Query_full.vcf.gz
tabix -p vcf Reference_full.vcf.gz
```
Now we use sample_full.map, and genetic_map still the same. --force-samples allows skiping undetected samples.

We run rfmix as follows:

```bash
./rfmix/rfmix -f Query_full.vcf.gz -r Reference_full.vcf.gz -m sample_full.map -g genetic_map_hg18_rfmix_modified.txt -o local15_PR_full --chromosome=chr15 --n-threads=20
```
We can change the number of --n-threads based on the available cores.

After the analysis we will obtain the same results as those in the pilot, lets see some results:

![4fbb4702-a00d-4da7-9f56-52cdceb087db](https://github.com/SilviaBuonaiuto/PANMIX/assets/55600771/befced81-5d41-4e36-9056-b77a8756078a)

IC95 is plotted, considering the number of samples. Because we have only two populations, we show ASW (1-ASW = CEU). 
We can analyze the density plot of the proportion:

![11230af9-b0b9-442d-9744-1fce593e00d9](https://github.com/SilviaBuonaiuto/PANMIX/assets/55600771/1353b1af-e71e-4fec-8dc0-768d19b18127)

This could serve us to detect unusual proportions of ancestral elements in the population, just a very preliminar analysis for developing hypothesis of natural selection.


# Running a full analysis (1 chromosome - five source superpopulations)
We performed some experiments in order to get a computational cost estimation. We analyzed the possibility of doing local ancestry with multiple source populations. Again we used Puerto Rico people from 1000gp as the query population. Also, chromosome 15 was analyzed.

One of the main limits was the RAM. With more than 7 source population we blow up a 64gb RAM with 32 cores (using just 5 cores) computer. In that case 543 reference (from source population) samples were used, and 134 query samples. 

We were able to compute local ancestry using the same number of reference samples, but with their superpopulation category (so only five superpopulations were available).

The computation time was 1.23 hours using 10 cores, RAM peaked at 43GB.

Here some results:

![Pop2](https://github.com/SilviaBuonaiuto/PANMIX/assets/55600771/1cab7605-897c-4c93-b5b3-84d6e8483119)


And if we want to analyze in more detail the variance in local ancestry:

![Pop](https://github.com/SilviaBuonaiuto/PANMIX/assets/55600771/2694da83-1505-4f01-b9eb-b217d8d2b095)


These results are consistent with those previously reported. Consider that CEU superpopulation is EU, and ASW superpopulation is AFR. 

# Running the COMPLETE analysis, all chromosomes with four source populations

For running this script you need being in a specific folder with the following scripts/files:
- rfmix: installed with executable
- sample_superpopulation.map: it contains four source populations from different continets consiering ror Europeans CEU and TSI populations, for Africans YRI MSL GWD populations, for Americans PEL, and for EAS the CHB.
- genetic_map_hg18_rfmix_modified.txt: downloaded with the link provided in this readme.

This script will produce 22 folders and download reference phased variant calling from 1000gp. Then it will extract reference samples.

```bash
#!/bin/bash

# Maximum number of processes to run in parallel
MAX_JOBS=22
jobs=0

for i in {1..22}; do
    (
        # Create a folder named 'chr' followed by the number
        mkdir "chr$i"

        # Change to the created directory
        cd "chr$i"

        # Download the genome files for each chromosome
        wget "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        wget "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

        # Rename the downloaded files
        mv "1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" "AllPop_chr$i.vcf.gz"
        mv "1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" "AllPop_chr$i.vcf.gz.tbi"

        # Use bcftools and tabix as before
        bcftools view -S ../reference_list.txt "AllPop_chr$i.vcf.gz" -Oz -o "Reference_chr$i.vcf.gz" --force-samples
        tabix -p vcf "Reference_chr$i.vcf.gz"

        # Return to the parent directory
        cd ..
    ) &

    # Control concurrency
    ((jobs=jobs+1))
    if [[ $jobs -ge $MAX_JOBS ]]; then
        wait -n
        ((jobs=jobs-1))
    fi
done

# Wait for all jobs to finish
wait
```

Once executed, query phased vcf should be splitted for different chromosomes, and each one should we moved to the specific folder, with the name: Query_chr$i.vcf.gz , being $i the number of the chromosome. Then we can run the analysis as follows:

```bash
#!/bin/bash

# Maximum number of processes to run in parallel
MAX_JOBS=22
jobs=0

for i in {1..22}; do
    (
        # Change to the directory corresponding to each chromosome
        cd "chr$i"

        # Run rfmix for local ancestry inference
        ./../rfmix -f "Query_chr$i.vcf.gz" -r "Reference_chr$i.vcf.gz" -m ../sample_superpopulation.map -g ../genetic_map_hg18_rfmix_modified.txt -o "local_chr$i" --chromosome="chr$i" --n-threads=20

        # Return to the parent directory
        cd ..
    ) &

    # Control concurrency
    ((jobs=jobs+1))
    if [[ $jobs -ge $MAX_JOBS ]]; then
        wait -n
        ((jobs=jobs-1))
    fi
done

# Wait for all jobs to finish
wait
```
