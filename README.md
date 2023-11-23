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








