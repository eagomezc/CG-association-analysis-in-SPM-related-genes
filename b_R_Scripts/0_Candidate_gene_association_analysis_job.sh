#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 5      # Request 10 core
#$ -l h_rt=72:0:0 # Request 20 hour runtime
#$ -l h_vmem=60G   # Request 200GB RAM
#$ -l highmem      # High Memory flag
#$ -t 1-26	# Number of lines in txt file	

# Create a list with all the chromosomes, genes and position we are interested:

chr=$(sed -n "${SGE_TASK_ID}p" chr.txt)
gene=$(sed -n "${SGE_TASK_ID}p" gene.txt)
start=$(sed -n "${SGE_TASK_ID}p" start.txt)
end=$(sed -n "${SGE_TASK_ID}p" end.txt)

# Create a for to interact and do the same for every chromosome:

module load plink

# Create a PED file only containing the samples we are interested in (Patients and Controls):

awk -F'\t' -v OFS="\t" 'NR==FNR{a[$1]; next} FNR==1 || $2 in a' ../../input/6_ped_manipulation/temp/For_RA_join.tab ../../input/6_ped_manipulation/ukb_chr${chr}.ped > ../../input/6_ped_manipulation/temp/ukb_chr${chr}_${gene}_RA_filter.ped

# Update the PED file so it will have a column with the different conditions (2 for cases and 1 for controls):

join -1 1 -2 1 -t $'\t' <(sort -k 1 ../../input/6_ped_manipulation/temp/For_RA_join.tab) <(sort -k 1 ../../input/6_ped_manipulation/temp/ukb_chr${chr}_${gene}_RA_filter.ped) | awk -F'\t' -v OFS="\t" '$7=$2' | cut -f2 --complement --output-delimiter=$'\t' > ../../output/6_ped_manipulation/ukb_chr${chr}_${gene}_RA.ped 

# Delete the PED file without the conditions:

rm ../../input/6_ped_manipulation/temp/ukb_chr${chr}_${gene}_RA_filter.ped

# Open the Quality Control folder:

cd ../../input/7_quality_controls/

# Creates soft link of the PED file that has only the correct samples with their conditions:

ln -s ../../output/6_ped_manipulation/ukb_chr${chr}_${gene}_RA.ped .
ln -s ../../output/5_from_bed_to_ped/ukb_chr${chr}.map ukb_chr${chr}_${gene}_RA.map

# Run PLINK to do some quality steps. First one is Missing Genotyping: 

plink --file ukb_chr${chr}_${gene}_RA --test-missing --out temp/ukb_chr${chr}_${gene}_RA

# Plink release a file that can be used to generate a list of SNPs that needs to be delete:

cd temp/

perl run-diffmiss-qc-name.pl ukb_chr${chr}_${gene}_RA

# Checks if there is duplicates or relatives in the samples. First it needs to be created a list with the samples id that can be obtained from the For_Condition_join.tab file:

awk '{print $1}' ../../../input/6_ped_manipulation/temp/For_RA_join.tab > list_id_RA.tab

# Remove_samples.py is a python script that generates a list of samples that needs to be delete: 

python remove_samples_RA.py

# Plink to run other quality controls: Remove relatives samples, failing SNPS, MAF, GENO, HWE. 

cd ..

# I used here a string threshold for HWE since we are working with autosomal chromosome. HWE in this case only considers controls.

plink --file ukb_chr${chr}_${gene}_RA --remove temp/exclude_samples_RA.tab --exclude temp/ukb_chr${chr}_${gene}_RA-fail-diffmiss-qc.txt -mind 0.1 --maf 0.01 --geno 0.2 --hwe 0.001 --recode --out ../../output/7_quality_controls/clean_ukb_chr${chr}_${gene}_RA

rm ../../output/6_ped_manipulation/ukb_chr${chr}_${gene}_RA.ped
rm ukb_chr${chr}_${gene}_RA.map

# Open the Statistical analysis folder:

cd ../../input/8_statistical_analysis/

# Creates soft link for the clean PED and MAP files: 

ln -s ../../output/7_quality_controls/clean_ukb_chr${chr}_${gene}_RA.ped .
ln -s ../../output/7_quality_controls/clean_ukb_chr${chr}_${gene}_RA.map .

# Creates a file with a column with all the samples surviving the QC step:

awk '{print $1}' clean_ukb_chr${chr}_${gene}_RA.ped > clean_samples_list_chr${chr}_${gene}_RA.tab

# Subset the standard covariates file so only contains the same samples as the PED file:

awk -F'\t' -v OFS="\t" 'NR==FNR{a[$1]; next} FNR==1 || $2 in a' clean_samples_list_chr${chr}_${gene}_RA.tab covariates_dummy_no_sex.tab > clean_cov_chr${chr}_${gene}_RA.tab

# Creates a directory that is going to keep all the results for that specific chromosome:

mkdir ../../output/8_statistical_analysis/RA_vs_Control/${chr}_all

mkdir ../../output/8_statistical_analysis/RA_vs_Control/${chr}_all/${gene}

# Run the association analysis and the models to identify relevant SNPs: 

# Adding "sex" in logistic specify sex as a covariate, so there is no need to add it in the covariates file. 
# Keep-pheno-on-missing-cov specify not to delete samples when there is no covariates info (useful for analysis with low sample number)
# Logistic and lineal when genotypic, hethom, dominant or recessive excludes automatically the non autosomal variants (meaning X chromosome variants).

plink --file clean_ukb_chr${chr}_${gene}_RA --logistic sex hide-covar --covar clean_cov_chr${chr}_${gene}_RA.tab keep-pheno-on-missing-cov --chr ${chr} --from-kb ${start} --to-kb ${end} --adjust --freq counts --out ../../output/8_statistical_analysis/RA_vs_Control/${chr}_all/${gene}/log_chr${chr}_${gene}_RA

plink --file clean_ukb_chr${chr}_${gene}_RA --logistic dominant sex hide-covar --covar clean_cov_chr${chr}_${gene}_RA.tab keep-pheno-on-missing-cov --chr ${chr} --from-kb ${start} --to-kb ${end} --adjust --freq counts --out ../../output/8_statistical_analysis/RA_vs_Control/${chr}_all/${gene}/log_chr${chr}_${gene}_RA_dom

plink --file clean_ukb_chr${chr}_${gene}_RA --logistic recessive sex hide-covar --covar clean_cov_chr${chr}_${gene}_RA.tab keep-pheno-on-missing-cov --chr ${chr} --from-kb ${start} --to-kb ${end} --adjust --freq counts --out ../../output/8_statistical_analysis/RA_vs_Control/${chr}_all/${gene}/log_chr${chr}_${gene}_RA_rec

# Remove in between ped and map files:

rm ../../output/7_quality_controls/clean_ukb_chr${chr}_${gene}_RA.ped
rm ../../output/7_quality_controls/clean_ukb_chr${chr}_${gene}_RA.map

# Remove in between files: 

rm clean_samples_list_chr${chr}_${gene}_RA.tab
rm clean_cov_chr${chr}_${gene}_RA.tab





