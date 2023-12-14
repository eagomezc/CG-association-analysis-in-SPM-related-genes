# ----------------- METAL ANALYSIS RESULTS  ---------------------------------#

# This scripts takes the results from CG analysis and the summary statistics
# from GWAS repositories to do a meta-analysis (combining the results from different
# studies to increase statistical power, validate results and confirm effect
# directionality). 

# The meta-analysis is done based on different algorithm: 

# Inverse-variance-weighted fixed effect - SCHEME STDERR:
  # Summarizes effect sizes from multiple independent studies by calculating the weighted 
  # mean of the effect sizes using the inverse variance of the individual studies as weights.
  # In other words, calculate the average of the effect having in consideration as weight for 
  # each individual value the variance. Meaning that more precise effects (less variance) have
  # more weight than less precise effects.
  # If the variances of the measurements are all equal, then the inverse-variance weighted 
  # average becomes the simple average.

# Sample-size-weighted - You need the number of samples for this one:
  # Constructs a new z-score by calculating a weighted sum of individual z-scores. It has been
  # known that the sample size of individual studies is a preferable weight for the method. 

# IMP:  IVW aims to obtain the best summary estimate (thus MLE and minimum variance), and SZ aims 
# to maximize the statistical power, without considering the summary estimator.

# Then it creates tables with the significant SNPs with their relevant information and plot the 
# results in a odd-ratio format plot. 

# ---> LIBRARY LOAD:

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

# To creates the QQplot of the models:

gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain((p-value))))
  log10Po <- expression(paste("Observed -log"[10], plain((p-value))))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 20, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, color = "red", size = 0.8) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

# To calculate the genetic inflation of the model. 
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

chr_color <- function(genes) {
 colors <- c("darkgoldenrod", "steelblue2", "violet", "olivedrab3",
             "palevioletred4", "plum2", "bisque4", "saddlebrown", 
             "orangered", "mediumorchid2", "lightskyblue", "mediumpurple", 
             "coral", "darkgoldenrod4", "lightskyblue4", "darkslateblue",
             "darkorchid1", "firebrick2", "lightgoldenrod3", "violetred1", 
             "orange2", "burlywood", "lightpink4", "chocolate2",
             "cornflowerblue", "darkorchid1", "gray49")  
 names(colors) <- c("GSTM4", "mGST3", "PTGS2", "LGR6", "mGST2", "LTC4S", "GPR37", "EPHX2", "ALOX5",
                "CMKLR1", "CYSLTR2", "GPR99", "GPR18", "LTB4R", "DPEP3", "DPEP2", "DPEP1", "ALOX15",
                "ALOX12", "ALOX15B", "GPR32", "FPR2", "GGT2", "GGT1", "CYSLTR1", "GPR101", "GPR173")
 col_sel <- colors[names(colors) %in% genes]
 col_sel
}

# ---> SET DIRECTORY: 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---> DATA MANIPULATION:

# Folder name with the results:

folders <- "METAL_RA_Results"
metal_names <- c("metal_RA_with_MegaGWAS_HE_SE_1", "metal_RA_with_MegaGWAS_HE_Ncases_1", "metal_RA_with_MegaGWAS_HE_SE_1")

# Open table with UK Biobank and validation results:

ukbiobank_results <- read.table(file = "METAL_RA_Results/all_SNPs_RA_vs_Control.tsv", 
                                header = TRUE,
                                sep = "\t",
                                stringsAsFactors = FALSE)

# Get only columns I care:
ukbiobank_results <- ukbiobank_results[, c(1:3, 7, 9, 10)]

Eur_all_final <- read.table(file = "METAL_RA_Results/Eur_all_final_Inc_chro_X_SPM.txt", 
                                header = TRUE,
                                sep = "\t",
                                stringsAsFactors = FALSE)

# Get only columns I care:
Eur_all_final <- Eur_all_final[, c(1:5)]

GWAS_RA <- read.table(file = "METAL_RA_Results/MegaGWAS_summary_European_SPM.txt", 
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE)

# Get only columns I care:
GWAS_RA <- Eur_all_final[, c(1:5)]

ukbiobank_results_rel <-  read.table(file = "METAL_RA_Results/rel_based_on_ngenes_SNPs_RA_vs_Control.tsv", 
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE)

# Get only columns I care:
ukbiobank_results_rel <- ukbiobank_results_rel[, c(1:3, 7, 9, 10)]

# ---> OPEN FILES:

for (z in 1: length(metal_names)) {
  
  meta_results <- read.table(file = paste(folders, "/", metal_names[[z]], ".tbl", sep = ""), 
                             header = TRUE,
                             sep = "\t",
                             stringsAsFactors = FALSE)
  
  meta_results$SNP <- meta_results$MarkerName # Get the SNPs columns as SNP.
  meta_results$MarkerName <- NULL
  
  # Here we only get the SNPs that still passed the p.value threshold (26 because # of genes that we study)
  # It can be down to p.value divided by number of SNPS identified in the UK biobank dataset maybe??????
  # Also, we want to identy the SNPs common in both databases which are the one with the calcaluated meta effect. 
  # When a SNP is in one but not in the other file, METAL just print the original values. 
  meta_rel_ngenes <- meta_results[meta_results$P.value <= 0.05/26 & !grepl("?", meta_results$Direction, fixed = TRUE), ]
  
  # Since we used the two METAL methods we need to identify which one has the Beta value to calculate the OR. 
  if("Effect" %in% colnames(meta_rel_ngenes)) {
  meta_rel_ngenes$OR_final <- exp(meta_rel_ngenes$Effect) } 
  
  # Here we got the info for every SNP by merging the meta result table with the ukbiobank results table.
  # Since this table only got SNPs from both analysis, the Uk biobank results is enough for this. 
  metal_rel_final <- merge(meta_rel_ngenes, ukbiobank_results, by = "SNP")
  metal_rel_final <- metal_rel_final[order(metal_rel_final$P.value), ] # Order for p.value. 
  
  # Save this table:
 # write.table(metal_rel_final, 
          #    file = paste("../../output/4_METAL/METAL_RA_Results/29_11_22_table_", metal_names[[z]], ".tsv", sep = ""),
           #   sep = "\t",
           #   quote = FALSE,
           #   row.names = FALSE)  
  
  # Get a table with only the significant SNPs that survive the Meta analysis from the original UK Biobank sig
  # results. 
  metal_rel_uk_final <- merge(meta_rel_ngenes, ukbiobank_results_rel, by = "SNP")
  metal_rel_uk_final <- metal_rel_uk_final[order(metal_rel_uk_final$P.value), ]
  
  # Save this table:
  #write.table(metal_rel_uk_final, 
             # file = paste("../../output/4_METAL/METAL_RA_Results/17_02_23_sig_table_", metal_names[[z]], ".tsv", sep = ""),
             # sep = "\t",
             # quote = FALSE,
             # row.names = FALSE) 
  
  # CREATE FIGURES:
  
  # Gene order:
  genes_order <- c("GSTM4", "mGST3", "PTGS2", "LGR6", "mGST2", "LTC4S", "GPR37", "EPHX2", "ALOX5",
                   "CMKLR1", "CYSLTR2", "GPR99", "GPR18", "LTB4R", "DPEP3", "DPEP2", "DPEP1", "ALOX15",
                   "ALOX12", "ALOX15B", "GPR32", "FPR2", "GGT2", "GGT1", "CYSLTR1", "GPR101", "GPR173")
  
  # Gene colors 
  color_genes <- chr_color(unique(metal_rel_final$ASSO_GENE))
  
  # So we want to have in the y-axis repeat label names, otherwise SNPs from the same gene will overlap in the 
  # same row. For that we just create a column with genes names and adding numbers so they will be different. 
  metal_rel_final$number <- c(1:nrow(metal_rel_final)) 
  metal_rel_final$label <- paste(metal_rel_final$ASSO_GENE, metal_rel_final$number, sep = "_")
  # Order based on genes and genome location:
  metal_rel_final$order <- match(metal_rel_final$ASSO_GENE, genes_order)
  metal_rel_final <- metal_rel_final[order(metal_rel_final$order), ]

  # Get figure y axis order: 
  metal_rel_final$label <- factor(metal_rel_final$label, levels = metal_rel_final$label)
  
  # Small loop to get repeat colors as well. The original function just save each color once. 
  color_y <- vector()
  for (c in 1:nrow(metal_rel_final)) { 
    color_add <- as.character(chr_color(metal_rel_final$ASSO_GENE[[c]]))
    color_y <- append(color_y, color_add)
  }
  
  # Get plots based on the type of table you have. 
  if("Effect" %in% colnames(metal_rel_final)) {
  
  beta <- 
    ggplot(metal_rel_final, aes(x=Effect, y=label)) +
    geom_point(aes(size=-log10(P.value), color=as.factor(ASSO_GENE))) + 
    geom_errorbarh(aes(xmin=Effect-StdErr, xmax=Effect+StdErr), position=position_dodge(.9)) + # All points
    scale_color_manual(values = color_genes) + 
    xlim(c(min(metal_rel_final$Effect)-0.05, max(metal_rel_final$Effect)+0.05)) + 
    scale_y_discrete(limits = rev, label = rev(metal_rel_final$ASSO_GENE)) +
    {if (nrow(metal_rel_final[!metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ]) > 0) geom_text(
      data = metal_rel_final[!metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ], 
              aes(max(metal_rel_final$Effect)+0.05, label, label = SNP), hjust = 0, size = 4, color = "black",
              position=position_dodge(0.5))} +
    {if (nrow(metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ]) > 0) geom_text(
      data = metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ], 
              aes(max(metal_rel_final$Effect)+0.05, label, label = SNP), color = "red", hjust = 0, 
              size = 4, position=position_dodge(0.5))} +
    geom_vline(xintercept = 0, color="black", linetype="dashed", size = 1) +
    guides(colour = "none") +
    labs(x = "Beta", y = "Genes", size = expression(-log[10](P))) +
    theme(
      axis.title = element_text(size = 20, colour = "black"),
      axis.text.x  =  element_text(size = 15, vjust = 0.5, colour = "black"), 
      axis.text.y  = element_text(size = 15, hjust = 1, colour = rev(color_y)),
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks = element_line(colour = "black", size = 0.5),  
      axis.ticks.length = unit(0.1, "cm"),
      legend.position="bottom",
      legend.text = element_text(size = 15, vjust = 0.5, colour = "black"), 
      legend.title = element_text(size = 15, vjust = 0.5, colour = "black"),
      legend.key = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray89", size = 0.5, linetype = 2),
      panel.background = element_blank())
  
  png(file = paste("../../output/4_METAL/METAL_RA_Results/17_02_23_Beta_", metal_names[[z]], ".png", sep = ""),
      width = 1200, height = 900)
  
  print(beta)
  
  dev.off()
  
  
  OR <- ggplot(metal_rel_final, aes(x=OR_final, y=label)) +
    geom_point(aes(size=-log10(P.value), color=as.factor(ASSO_GENE))) + 
    geom_errorbarh(aes(xmin=OR_final-StdErr, xmax=OR_final+StdErr), position=position_dodge(.9)) + # All points
    scale_color_manual(values = color_genes) + 
    xlim(c(min(metal_rel_final$OR_final)-0.05, max(metal_rel_final$OR_final)+0.05)) + 
    scale_y_discrete(limits = rev, label = rev(metal_rel_final$ASSO_GENE)) +
    {if (nrow(metal_rel_final[!metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ]) > 0) geom_text(
      data = metal_rel_final[!metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ], 
              aes(max(metal_rel_final$OR_final)+0.05, label, label = SNP), hjust = 0, size = 4, color = "black",
              position=position_dodge(0.5))} +
    {if (nrow(metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ]) > 0) geom_text(
      data = metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ], 
              aes(max(metal_rel_final$OR_final)+0.05, label, label = SNP), color = "red", hjust = 0, 
              size = 4, position=position_dodge(0.5))} +
    geom_vline(xintercept = 1, color="black", linetype="dashed", size = 1) +
    guides(colour = "none") +
    labs(x = "Odd ratios", y = "Genes", size = expression(-log[10](P))) +
    theme(
      axis.title = element_text(size = 20, colour = "black"),
      axis.text.x  =  element_text(size = 15, vjust = 0.5, colour = "black"), 
      axis.text.y  = element_text(size = 15, hjust = 1, colour = rev(color_y)),
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.ticks = element_line(colour = "black", size = 0.5),  
      axis.ticks.length = unit(0.1, "cm"),
      legend.position="bottom",
      legend.text = element_text(size = 15, vjust = 0.5, colour = "black"), 
      legend.title = element_text(size = 15, vjust = 0.5, colour = "black"),
      legend.key = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray89", size = 0.5, linetype = 2),
      panel.background = element_blank())
  
  png(file = paste("../../output/4_METAL/METAL_RA_Results/17_02_23_OR_", metal_names[[z]], ".png", sep = ""),
      width = 1200, height = 900)
  
  print(OR)
  
  dev.off()  } else if("Zscore" %in% colnames(metal_rel_final)) {
        
      Zscore <- 
          ggplot(metal_rel_final, aes(x=Zscore, y=label)) +
          geom_point(aes(size=-log10(P.value), color=as.factor(ASSO_GENE))) + 
          scale_color_manual(values = color_genes) + 
          xlim(c(min(metal_rel_final$Zscore)-2, max(metal_rel_final$Zscore)+2)) + 
          scale_y_discrete(limits = rev, label = rev(metal_rel_final$ASSO_GENE)) +
         {if(nrow(metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ]) > 0) 
          geom_text(data = metal_rel_final[metal_rel_final$SNP %in% metal_rel_uk_final$SNP, ], 
                    aes(max(metal_rel_final$Zscore)+1.5, label, label = SNP), color = "red", hjust = 0, 
                    size = 4, position=position_dodge(0.5)) } + 
          geom_vline(xintercept = 0, color="black", linetype="dashed", size = 1) +
          guides(colour = "none") +
          labs(x = "Zscore", y = "Genes", size = expression(-log[10](P))) +
          theme(
            axis.title = element_text(size = 20, colour = "black"),
            axis.text.x  =  element_text(size = 15, vjust = 0.5, colour = "black"), 
            axis.text.y  = element_text(size = 15, hjust = 1, colour = rev(color_y)),
            axis.line = element_line(colour = 'black', size = 0.5), 
            axis.ticks = element_line(colour = "black", size = 0.5),  
            axis.ticks.length = unit(0.1, "cm"),
            legend.position="bottom",
            legend.text = element_text(size = 15, vjust = 0.5, colour = "black"), 
            legend.title = element_text(size = 15, vjust = 0.5, colour = "black"),
            legend.key = element_rect(fill = "white"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(color = "gray89", size = 0.5, linetype = 2),
            panel.background = element_blank())
        
        
        png(file = paste("../../output/4_METAL/METAL_RA_Results/17_02_23_Zscore_", metal_names[[z]], ".png", sep = ""),
            width = 1200, height = 900)
        
        print(Zscore)
        
        dev.off() 
  
      } else {next}
}
