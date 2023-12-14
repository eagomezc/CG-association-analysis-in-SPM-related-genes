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

## [a_Data](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/a_Data)

This folder contains, separated by subfolders, the different file formats that has to be used to run the different R scripts. 

**NOTE:** Genotype and phenotype information from UK Biobank is not publicly available. 

The subfolders are:

**1_Visualization_plots_(R_script)**: Contains a folder with summary statistics for all the SPM-related genes. In addition, it contains tables with gene's information (chromosome and gene name) and available information of the genes in UK Biobank. 

**2_Meta-analysis_(METAL)**: Contains summary statistics files for meta-analysis in addition of the configuration file ([metal_RA_sig_with_MegaGWAS_sig_HE_SE.txt](https://github.com/eagomezc/CG-association-analysis-in-SPM-related-genes/blob/main/a_Data/2_Meta-analysis_(METAL)/metal_RA_sig_with_MegaGWAS_sig_HE_SE.txt)) to run METAL.

**3_METAL_visualization_results**: Contains tab-delimited results from META analysis. 

More details about the format of this files can be seen in the comments of each script. 

## [b_R_Scripts](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/b_R_Scripts)

This folder contains the scripts used to create the support vector machine and random forest prediction models, in addition to the script used to run differential gene expression analysis of interested genes in the lipid mediator pathways. 

The scripts are: 

**1_machine_learning_(All_methodologies).R**: Using a training dataset, this script creates the machine learning models (Bayesian, Elastic net, SVM and random forest) used to predict the response to DMARD treatment in rheumatoid arthritis patient. The script works with the packages **classyfire, randomForest, glmnet, caret, arm** that use different machine learning methodologies and bootstrapping for the model creation. 

**2_randomForest_(RF_models).R**: Using a training dataset, this script creates the machine learning models used to predict the response to DMARD treatment in rheumatoid arthritis patient. The script works with the package **randomForest** that uses random forests and bootstrapping for the model creation. Besides that, estimate the **importance** of each lipid mediator in the improvement of the model's accuracy. Finally, it also uses the test cohort to evaluate the models and estimate the area under the receiver operating characteristic curves (AUC). 

**3_DGE_analysis_(Edge_R).R**: Using RNA-seq raw read counts, this scripts performs differential gene expression analysis using the package **Edge R**, that uses the quasi-likelihood method to identify differences in the expression levels of specific genes between the DMARD responder and Non Responder rheumatoid arthritis patients. It also creates violin plots as a way to visualize the different gene expression levels.

More details of how the scripts works can be seen in the comments of each script. 

## [c_Expected_Output](https://github.com/eagomezc/Machine-Learning-and-RA-treatment/tree/main/c_Expected_Output)

This folder contains, separated by subfolders, the different expected outputs that can be obtain after running the R scripts. Each subfolder has the name of the specific script that generates it, in addition to the number of the script, to make more clear what file is the result of what script. At the moment to download this repository in a local computer, it's important to remember that all the **output pathways in the scripts has to be changed**.

The subfolders are:

**1_machine_learning_(All_methodologies).R**: The expected results from this script are a tab-delimited file containing a table with the model's names, the machine learning strategy used, their accuracy percentages, sensitivity, specificity and confusion table; a figure of all the models with their accuracy score, the tunning parameters figure and the different models saved as an R object that can be used in the future.  

**2_randomForest_(RF_models)**: The expected results from this script are a tab-delimited file containing a table with the model's names, their accuracy percentages and their AUC values after evaluation with the test cohort; the different models saved as an R object that can be used in the future; and pdf files that contains plots associated with the performance of the models and the importance of each lipid mediator in the construction of the models. 

**3_DGE_analysis_(Edge_R)**: The expected results from this script is a tab-delimited file containing a table with the gene's names, their log(FC), log(CPM), F value, p value and adjust p value (FDR). In addition, is expected to generate a pdf file with violin plots of ALOX-related enzymes.

More details about how this files are generated can be seen in the comments of each script. 

# Publication:

Part of the results from this section of my thesis are described in the following paper: 

[Gomez, E.A., Colas, R.A., Souza, P.R., Hands, R., Lewis, M.J., Bessant, C., Pitzalis, C., Dalli, J., 2020. Blood pro-resolving mediators are linked with synovial pathology and are predictive of DMARD responsiveness in rheumatoid arthritis. Nat Commun 11, 5420.](https://www.nature.com/articles/s41467-020-19176-z) 
 
 





