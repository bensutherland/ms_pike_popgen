# Northern pike population genetics study
A repository to accompany the analysis of the northern pike popgen study.     

Citation:     
```
Johnson, HA, Rondeau, EB, Minkley, DR, Leong, JS, Whitehead, J, Despins, CA, et al. (2020). Population genomics of North American northern pike: variation and sex-specific signals from a chromosome-level, long read genome assembly. bioRxiv, 2020.2006.2018.157701. doi:10.1101/2020.06.18.157701

```

#### Requirements ####
- samtools     
- bedtools     
- bcftools
- vcftools
- Eric Normandeau's scripts repository

## 00. Data preparation
Download the following and put in `02_input_data`:       
- _Esox lucius_ reference genome [GCA_004634155.1_Eluc_v4_genomic.fna](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004634155.1://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004634155.1/)       
- project VCF file from [figshare](#TODO)     
 

## 01. Characterize overall variant counts across chromosomes
Prepare and characterize the chromosomes:     
```
# Change directory
cd 02_input_data

# Trim scaffold names to match formats observed in the VCF e.g., LG03    
sed 's/.*LG/>LG/g' GCA_004634155.1_Eluc_v4_genomic.fna | sed 's/,.*//g' > GCA_004634155.1_Eluc_v4_genomic_renamed.fna    

# Unwrap fasta
fasta_unwrap.py GCA_004634155.1_Eluc_v4_genomic_renamed.fna GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap.fna

# Keep chromosomes only
grep -E -A1 '^>LG' GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap.fna | grep -vE '^--$' - > GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap_chr_only.fna  

# Determine lengths of chromosomes
`fasta_lengths.py GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap_chr_only.fna`    

# Index the reference
samtools faidx ./GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap_chr_only.fna 

```

Determine windows positions per chromosome and extract variants:      
```
# Prepare a windows file      
bedtools makewindows -g ./GCA_004634155.1_Eluc_v4_genomic_renamed_unwrap_chr_only.fna.fai -w 1000000 > windows.bed

# Extract variants per window from each chromosome
bedtools coverage -a windows.bed -b ./Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode.vcf -counts > coverage.txt    
```

Use the following script interactively in R to analyze the variants per window and generate plots:      
`01_scripts/plotting_chr_poly_levels.R`    


## 02. Identify fixed differences across and within population clusters ##
Use the following script interactively in R to analyze fixed genotypic differences across and within the major groupings in the dataset:     
`01_scripts/fixed_diffs_across_NACD.R` .    


## 03. Per sample genotypes
Produce a stacked barplot of genotypes per sample, including REF/REF, REF/ALT, ALT/ALT, missing using the following script:      
`01_scripts/per_sample_genotypes.R`      


## 04. Per population statistics
Variant sites with more than two alleles are causing issues, so first, remove these:     
`bcftools view --max-alleles 2 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode.vcf > 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf`        
This file should now have 1,117,361 variants remaining.     

Use the RScript `01_scripts/per_pop_stats.R` interactively to calculate per population average HOBS and FIS, and the number of private alleles per population.      


### Tajima's D ###
Using the VCF above with only biallelic variants, determine the names of the samples in the VCF:     
`bcftools query --list-samples 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf > 02_input_data/samples.txt`     

Create per population sample names files:     
```
grep '\-CR' 02_input_data/samples.txt > 02_input_data/samples_CHT.txt
grep '\-YR' 02_input_data/samples.txt > 02_input_data/samples_HOO.txt
grep '\-PL' 02_input_data/samples.txt > 02_input_data/samples_PAL.txt
grep '\-Mb' 02_input_data/samples.txt > 02_input_data/samples_WHI.txt
grep '\-S' 02_input_data/samples.txt > 02_input_data/samples_SLA.txt
grep '\-NJ' 02_input_data/samples.txt > 02_input_data/samples_HCK.txt
```

Create population-specific VCFs:     
```
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_CHT.txt --recode --recode-INFO-all --out 02_input_data/Eluc_CHT.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_HOO.txt --recode --recode-INFO-all --out 02_input_data/Eluc_HOO.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_PAL.txt --recode --recode-INFO-all --out 02_input_data/Eluc_PAL.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_WHI.txt --recode --recode-INFO-all --out 02_input_data/Eluc_WHI.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_SLA.txt --recode --recode-INFO-all --out 02_input_data/Eluc_SLA.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_HCK.txt --recode --recode-INFO-all --out 02_input_data/Eluc_HCK.vcf
vcftools --vcf 02_input_data/Eluc.variants.GATK.iteration.2.b.Full.SNP.GATK_HF_removed.minQ20.mmdp10.mxmdp60.mmc10.mac1.homsumfilt.new.samp.names.recode_biallelic_only.vcf --keep 02_input_data/samples_greater_NA.txt --recode --recode-INFO-all --out 02_input_data/Eluc_GNA.vcf

```

Run Tajima's D calculation on the population-specific VCF:     
```
vcftools --vcf 02_input_data/Eluc_CHT.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_CHT
vcftools --vcf 02_input_data/Eluc_HOO.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_HOO
vcftools --vcf 02_input_data/Eluc_PAL.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_PAL
vcftools --vcf 02_input_data/Eluc_WHI.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_WHI
vcftools --vcf 02_input_data/Eluc_SLA.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_SLA
vcftools --vcf 02_input_data/Eluc_HCK.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_HCK
vcftools --vcf 02_input_data/Eluc_GNA.vcf.recode.vcf --TajimaD 10000 --out 03_results/Eluc_GNA

```
note: output will be, for example, `Eluc_CHT.Tajima.D`.    

Run Tajima's D calculation on the population-specific VCF:     
```
vcftools --vcf 02_input_data/Eluc_CHT.vcf.recode.vcf --site-pi --out 03_results/Eluc_CHT
vcftools --vcf 02_input_data/Eluc_HOO.vcf.recode.vcf --site-pi --out 03_results/Eluc_HOO
vcftools --vcf 02_input_data/Eluc_PAL.vcf.recode.vcf --site-pi --out 03_results/Eluc_PAL
vcftools --vcf 02_input_data/Eluc_WHI.vcf.recode.vcf --site-pi --out 03_results/Eluc_WHI
vcftools --vcf 02_input_data/Eluc_SLA.vcf.recode.vcf --site-pi --out 03_results/Eluc_SLA
vcftools --vcf 02_input_data/Eluc_HCK.vcf.recode.vcf --site-pi --out 03_results/Eluc_HCK
vcftools --vcf 02_input_data/Eluc_GNA.vcf.recode.vcf --site-pi --out 03_results/Eluc_GNA

```
Use Rscript `X.R` to find means for each of the Tajima's D files.    

