# Genetic variants and regulation of specialized pro-resolving mediator in rheumatoid arthritis

# Overview: 

This repository contains the **Bash and R scripts** used to do a candidate gene association analysis with the purpose of identifying genetic variants in SPM-related genes associated with rheumatoid arthritis. Genotype and phenotype data from participant was obtained from the [**UK Biobank dataset**](https://www.ukbiobank.ac.uk/) and association analysis were performed using [**Plink 1.9 tools**](https://www.cog-genomics.org/plink/).

In addition to this, I also used a couple of **R scripts** for summary statistics visualization and [**METAL software**](https://genome.sph.umich.edu/wiki/METAL_Documentation) for meta-analysis. 

**NOTE:** **FUMA GWAS web-based application** (more information [here](https://fuma.ctglab.nl/)) (version 1.6.0) was used for fine mapping and functional annotation of the candidate genetic variants (biological consequences). 

# System Requirements: 

## Hardware requirements: 

All the scripts and software used for the **PhD Thesis** were run in a standard computer (RAM: 8GB, CP$: 4 cores, 3.60 GHZ/core). 

Big data analysis (including UK Biobank data manipulation and association analysis) were performed in the high-performance computing cluster of QMUL (more information [here](https://docs.hpc.qmul.ac.uk/)). Running time in the cluster depends of memory availability. 

## System requirements:

All the scripts were created and used on **Windows 10**:

**R version**: 4.0.4 

**R Studio version**: 1.1.524

The scripts should be compatible with Mac and Linux operating systems. 

For installing R and R Studio, follows the installation instructions [here](https://www.stats.bris.ac.uk/R/) and [here](https://www.rstudio.com/products/rstudio/download/). With typical installation times on a normal computer not exceeding 3h.

Plink 1.9 tools are already installed in high-performance computing cluster and it was called using the function **module load plink/1.9-170906**.

For installing METAL, follows the installation instructions [here](https://csg.sph.umich.edu/abecasis/metal/download/).  

## Required R packages (libraries): 

The requiered packates to run all scripts should be installed automatically when you run each scritp; however, in case you need to install them manually, you can do it as follow:

```
# Packages ggplot2, ggrepel, dplyr, gggenes:
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('gggenes')) install.packages('gggenes'); library('gggenes')

```
# Content: 

The repository contains three folders: 

## [a_Data](https://github.com/eagomezc/CG-association-analysis-in-SPM-related-genes/tree/main/a_Data)

This folder contains, separated by subfolders, the different file formats that has to be used to run the different R scripts. 

**NOTE:** Genotype and phenotype information from UK Biobank is not publicly available. 

The subfolders are:

**1_Visualization_plots_(R_script)**: Contains a folder with summary statistics for all the SPM-related genes. In addition, it contains tables with gene's information (chromosome and gene name) and available information of the genes in UK Biobank. 

**2_Meta-analysis_(METAL)**: Contains summary statistics files for meta-analysis in addition of the configuration file ([metal_RA_sig_with_MegaGWAS_sig_HE_SE.txt](https://github.com/eagomezc/CG-association-analysis-in-SPM-related-genes/blob/main/a_Data/2_Meta-analysis_(METAL)/metal_RA_sig_with_MegaGWAS_sig_HE_SE.txt)) to run METAL.

**3_METAL_visualization_results**: Contains tab-delimited results from META analysis. 

More details about the format of this files can be seen in the comments of each script. 

## [b_R_Scripts](https://github.com/eagomezc/CG-association-analysis-in-SPM-related-genes/tree/main/b_R_Scripts)

This folder contains the scripts used to performed and analysed the candidate gene association analysis and meta-analysis. 

The scripts are: 

**0_Candidate_gene_association_analysis_job.sh**: Using the genotype and phenotype information of cases and controls from UK Biobank, this script perfomed association analysis of selected SPM-related genes. The script **filter the cases and the controls**, performed **quality controls** steps and finally run the association analysis based on **additive, recessive and dominant genetic models**. 

**1_Visualization_plots_(R_script).R**: Takes the summary statistics obtained from the previous script and generate a manhattan plot for visualization of genetic variants associated with rheumatoid arthritis.

**NOTE**: Meta-analysis using METAL are run using the terminal. To run METAL I used the following command: **../metal config.txt**. 

**3_METAL_visualization_results_(R_script).R**: Taking the meta-analysis results from the previous analysis, this script generates a forest plot highlighting the effect size of the candidate genetic variant.  

More details of how the scripts works can be seen in the comments of each script. 

## [c_Expected_Output](https://github.com/eagomezc/CG-association-analysis-in-SPM-related-genes/tree/main/c_Expected_Output)

This folder contains, separated by subfolders, the different expected outputs that can be obtain after running the R scripts. Each subfolder has the name of the specific script that generates it, in addition to the number of the script, to make more clear what file is the result of what script. At the moment to download this repository in a local computer, it's important to remember that all the **output pathways in the scripts has to be changed**.

The subfolders are:

**1_Visualization_plots_(R_script)**: The expected results from this script are candidate gene summarize statistics and manhattan plots for visualization. 

**3_METAL_visualization_results_(R_script)**: The expected results from this script are a tab-delimited file containing meta-analysis results including genetic variants' information (location, p-values, beta scores, odd-ratios, etc.) and Forest plot summarizing the results.  

More details about how this files are generated can be seen in the comments of each script. 

# Publication:

Part of the results from this section of my thesis are described in a puplication that is been peer-reviewed. 





