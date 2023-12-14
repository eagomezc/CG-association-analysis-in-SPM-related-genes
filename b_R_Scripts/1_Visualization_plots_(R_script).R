# ----------------- CG VISUALIZATION ANALYSIS ---------------------------------#

# This scripts takes the results from plink (CG analysis) and put together the
# results for all the genes, so FDR can be made to all of them at once.

# It creates Manhattan plots, calculates the genomic inflation factor, separates
# the SNPS of interest by gene and location within the gene, generates a list of
# candidate SNPs so it can be applied LD using the UK Biobank dataset and finally
# corroborates the SNPs with Gtex database so we have information about cis-eqtl.

# ---> LIBRARY LOAD:

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('gggenes')) install.packages('gggenes'); library('gggenes')

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

folders <- c("RA_vs_Control_NOiC")
condition <- "RA"

# ---> OPEN FILES:

# Upload file with chromosomes and genes names of interest: 

chr_genes <- read.table(file = "genes_CG.tsv", 
                        sep = "\t",
                        header = TRUE,
                        stringsAsFactors = FALSE)

ngenes <- nrow(chr_genes)


# File with all the SNPs information from Ukbiobank:

ukbiobank_snps <- read.table(file = "UK_biobank_full_SNP.tsv", 
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE)

# Validation results:

# This file comes from the comparison between the SNPs found in UK Biobank and those in the 
# validation cohorts. 

if (file.exists(paste(condition, "_final_Uk_and_validation_results_SNPs.tsv", sep = "")) == TRUE) {

val_results <- read.table(file = paste(condition, "_final_Uk_and_validation_results_SNPs.tsv", sep = ""), 
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)
}

# ---> DATA MANIPULATION:

chr <- c("1", "4", "5", "7", "8", "10", "12", "13", "14", "16", "17", "19", "22", "X")
asso <- c( "", "_dom", "_gen", "_hethom", "_rec")

for (z in 1: length(folders)) {

for (a in 1:length(asso)) {
  
  # Creates the tables with all the SNPs and P values for each association:
  
  all_table <- data.frame(CHR = 1,
                          SNP = "del",
                          BP = 1,
                          A1 = "del",
                          TEST = "del",
                          NMISS = 1,
                          OR = 1,
                          STAT = 1,
                          P = 1, 
                          ASSO_GENE = "gene")

  for (i in 1:length(chr)) {
  
  # Open the files for each gene:
  
  genes_in <- chr_genes[chr_genes$Chromosome == chr[[i]], ]
    
    for (k in 1:nrow(genes_in)) {
      
      # Check file exists:
      
      if (file.exists(paste(folders[[z]], "/", chr[[i]], "_all/", 
                            genes_in$Gene[k], "/log_chr", chr[[i]],
                            "_", genes_in$Gene[k], "_", condition, 
                            asso[[a]], ".assoc.logistic", sep = "")) == TRUE) {

      # Logistic association model 1:
      
      log_model <- read.table(file = paste(folders[[z]], "/", chr[[i]], "_all/", 
                                           genes_in$Gene[k], "/log_chr", chr[[i]],
                                           "_", genes_in$Gene[k], "_", condition, 
                                           asso[[a]], ".assoc.logistic", sep = ""), 
                              header = TRUE,
                              stringsAsFactors = FALSE)
      
      log_model$ASSO_GENE <- genes_in$Gene[k]
      
      all_table <- rbind(all_table, log_model) } else {next}

    }
    
  }
  
  all_table <- all_table[-1, ] # Remove toy first row.
  
  # Remove all the rows incomplete (meaning no P value or OR value, etc.)
  all_table <- na.omit(all_table) 
  
  # Order the table based on SNP and P value so when we delete duplicates, we 
  # will delete the duplicate with the higher P value. 
  all_table <- all_table[order(all_table$SNP, all_table$P), ]
  
  all_table <- all_table[!duplicated(all_table$SNP), ] # Duplication deletion.
  
  # Calculates the multiple comparison correction: 
  all_table$BONF <- p.adjust(all_table$P, "bonferroni")
  all_table$HOLM <- p.adjust(all_table$P, "holm")
  all_table$BH <- p.adjust(all_table$P, "BH")
  all_table$HOMMEL <- p.adjust(all_table$P, "hommel")
  all_table$FDR <- p.adjust(all_table$P, "fdr")
  
  # Get the mean of all the FDR corrections: 
  all_table$CON_ADJ <- rowSums(all_table[, c(11:15)])/5
  
  # Merge the SNP table with the info from UK Biobank:
  
  all_table_final <- merge(all_table, ukbiobank_snps, by = "SNP", all.x = TRUE)
  
  # Order the table by Chromosome and position. Easy for the graph making.
  all_table_final <- all_table_final[order(all_table_final$CHR, all_table_final$BP), ]
  all_table_final <- all_table_final[!duplicated(all_table_final$SNP), ]
  
  # Save the table with the name of the related model

  all_table_final_rel <- all_table_final[all_table_final$P <= 0.05, ] # Only significant SNPs
  all_table_final_rel <- all_table_final_rel[order(all_table_final_rel$CON_ADJ), ]
  
  all_table_final_rel_ngenes <- all_table_final_rel[all_table_final_rel$P <= (0.05/ngenes), ]
  
  if (nrow(all_table_final_rel) > 0) {
    
    # One table safe as a results:
    
    write.table(all_table_final_rel, 
              file = paste("../../output/1_SNPs_identification/", folders[[z]], "/", "relevant_SNPs", 
                           asso[[a]], "_", folders[[z]], ".tsv", sep = ""),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE) 
    
    write.table(all_table_final, 
                file = paste("../../output/1_SNPs_identification/", folders[[z]], "/", "all_SNPs", 
                             asso[[a]], "_", folders[[z]], ".tsv", sep = ""),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE) 
    
    write.table(all_table_final_rel_ngenes, 
                file = paste("../../output/1_SNPs_identification/", folders[[z]], "/", "rel_based_on_ngenes_SNPs", 
                             asso[[a]], "_", folders[[z]], ".tsv", sep = ""),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE) 
    
  }
  
  # ---> FIGURES AND DATA SEPARATION:
  
  # Get the tables and table format to make the Manhattan plot for every gene:
  
  man_table <- all_table_final 
  
  # Creates a new column with a sequence from 1 to nrow(man_table) to represent
  # relative location of the SNPs in the genome. 
  
  man_table$BPcum <- 1
  first <- 3
  
  for (m in unique(man_table$CHR)) {
    
   man_table[man_table$CHR == m, ]$BPcum <- seq(first, length.out = nrow(man_table[man_table$CHR == m, ]), by = 1)
   first <- max(man_table$BPcum) + 10
    
  }
  
  # Get a table with the SNPs that are associated with the gene based on the info
  # from UK Biobank. 
  
  table_genes <- subset(man_table, grepl(paste(chr_genes$Gene, collapse= "|"), 
                                         UK_ASSO_GENE))
  
  # Creates a new table with the information of where to locate the label for the
  # genes in the Manhattan plot. It will be the median of the sequence fraction
  # associated with that gene.
  
  label_pos <- man_table %>% 
    
    # Here you group by genes (a column with genes names) and then summarize creates
    # a new column with the desireable analysis. In this case the median between
    # the maximum and miminimum BPcum. 
    group_by(ASSO_GENE) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2)
  
  # Get the first gene from each chromosome to avoid overlapping:
 # first_gene <- chr_genes[!duplicated(chr_genes$Chromosome), ]
  
  # Delete the first gene that it's not problematic:
 # first_gene <- first_gene[-1, ]
  
  # Add plus 10 to the first gene position of each Chromosome to avoid overlapping
 # for (e in 1:nrow(label_pos)) {
 #  if (label_pos$ASSO_GENE[e] %in% first_gene$Gene) {
 #     label_pos$center[e] <- label_pos$center[e] + 0.8 }
 # }
    
  # Established the colors for the different genes:
  color_genes <- chr_color(unique(man_table$ASSO_GENE))
  
  # Established Y limit based on the lowest p value and BP cum.
  max_p_value <- -log10(min(man_table$P))
  
  # Established Y intersection where the highest p value pass multiple comparison
  # corrections: 
  
  if (nrow(man_table[man_table$BH <= 0.05, ]) > 0) { 
    fdr_intercept <- -log10(max(man_table[man_table$BH <= 0.05, ]$P)) - 0.1 } else {
      fdr_intercept <- -1 }
  
  # Creates purple table with common SNPs between validation and UK:
  
  if (file.exists(paste(condition, "_final_Uk_and_validation_results_SNPs.tsv", sep = "")) == TRUE) {
    
    man_table_purple <- merge(val_results, man_table, by = c("SNP", "CHR", "BP")) 
    
    # Creates the Manhattan plot:
    
    manhatan_plot <- ggplot(man_table, aes(x=BPcum, y=-log10(P))) +
      geom_point(aes(color=as.factor(ASSO_GENE)), size=3) + # All points
      geom_point(data = table_genes, aes(color=as.factor(ASSO_GENE)), alpha=1, shape=11 ,size=3) + # Only SNPs associated with genes.
      scale_color_manual(values = color_genes) +
      # Location of the labels for every gene:
      scale_x_continuous(label = label_pos$ASSO_GENE, breaks= label_pos$center) +
      # Add the panels for chromosome:
      facet_grid(.~ CHR, space = 'free_x', scales = 'free_x', switch = 'x') +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Chromosome") +
      coord_cartesian(ylim = c(0, max_p_value+0.3)) + #Y limit
      # Horizontal line representing significant SNPs:
      geom_hline(yintercept = fdr_intercept, color="red", linetype="dashed") +
      geom_hline(yintercept = -log10(0.05/ngenes), color="blue", linetype="dashed") +
      # Add SNP name to significant SNPs: 
      geom_text_repel(data = man_table[man_table$P <= (0.05/ngenes) & man_table$BH > 0.05, ], position = "identity",
                      label = man_table[man_table$P <= (0.05/ngenes) & man_table$BH > 0.05, ]$SNP, size = 4, color="blue",  
                      direction="x", angle=90, hjust=1, vjust=0.35) +
      geom_text_repel(data = man_table[man_table$BH <= 0.05, ], position = "identity", 
                      label = man_table[man_table$BH <= 0.05, ]$SNP, size = 4, color="red",  
                      direction="x", angle=90, hjust=1, vjust=0.35) +
      geom_text_repel(data = man_table_purple, position = "identity", 
                      label = man_table_purple$SNP, size = 4, color="purple",  
                      direction="x", angle=90, hjust=1, vjust=0.35, max.overlaps = Inf) +
      labs(caption = "Red line = Adj pvalue (BH) < 0.05 \n Blue line = pvalue < 0.05/number of study genes \n
         Purple SNPs = Match with Val cohort \n Stars = SNP in gene", color = "black", size = 3) +
      theme(
        axis.title = element_text(size = 15, colour = "black"),
        axis.text.x  =  element_text(size = 12, vjust = 0.5, angle =90, colour = "black"), 
        axis.text.y  = element_text(size = 12, hjust = 1, colour = "black"),
        axis.line = element_line(colour = 'black', size = 0.5), 
        axis.ticks = element_line(colour = "black", size = 0.5),  
        axis.ticks.length = unit(0.1, "cm"),
        legend.position="none",
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "ghostwhite"),
        strip.placement = 'outside',
        plot.caption = element_text(color = "Black", size = 10))
    
    
    png(file = paste("../../output/1_SNPs_identification/", folders[[z]], "/mp",  asso[[a]], ".png", sep = ""),
        width = 1200, height = 900)
    
    print(manhatan_plot)
    
    dev.off()
    
    
    } else {
    
      # Creates the Manhattan plot:
      
      manhatan_plot <- ggplot(man_table, aes(x=BPcum, y=-log10(P))) +
        geom_point(aes(color=as.factor(ASSO_GENE)), size=3) + # All points
        geom_point(data = table_genes, aes(color=as.factor(ASSO_GENE)), alpha=1, shape=11 ,size=3) + # Only SNPs associated with genes.
        scale_color_manual(values = color_genes) +
        # Location of the labels for every gene:
        scale_x_continuous(label = label_pos$ASSO_GENE, breaks= label_pos$center) +
        # Add the panels for chromosome:
        facet_grid(.~ CHR, space = 'free_x', scales = 'free_x', switch = 'x') +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Chromosome") +
        coord_cartesian(ylim = c(0, max_p_value+0.3)) + #Y limit
        # Horizontal line representing significant SNPs:
        geom_hline(yintercept = fdr_intercept, color="red", linetype="dashed") +
        geom_hline(yintercept = -log10(0.05/ngenes), color="blue", linetype="dashed") +
        # Add SNP name to significant SNPs: 
        geom_text_repel(data = man_table[man_table$P <= (0.05/ngenes) & man_table$BH > 0.05, ], position = "identity",
                        label = man_table[man_table$P <= (0.05/ngenes) & man_table$BH > 0.05, ]$SNP, size = 4, color="blue",  
                        direction="x", angle=90, hjust=1, vjust=0.35) +
        geom_text_repel(data = man_table[man_table$BH <= 0.05, ], position = "identity", 
                        label = man_table[man_table$BH <= 0.05, ]$SNP, size = 4, color="red",  
                        direction="x", angle=90, hjust=1, vjust=0.35) +
        # geom_text_repel(data = man_table_purple, position = "identity", 
        #                label = man_table_purple$SNP, size = 4, color="purple",  
        #                direction="x", angle=90, hjust=1, vjust=0.35, max.overlaps = Inf) +
        labs(caption = "Red line = Adj pvalue (BH) < 0.05 \n Blue line = pvalue < 0.05/number of study genes \n 
             Stars = SNP in gene", color = "black", size = 3) +
        theme(
          axis.title = element_text(size = 15, colour = "black"),
          axis.text.x  =  element_text(size = 12, vjust = 0.5, angle =90, colour = "black"), 
          axis.text.y  = element_text(size = 12, hjust = 1, colour = "black"),
          axis.line = element_line(colour = 'black', size = 0.5), 
          axis.ticks = element_line(colour = "black", size = 0.5),  
          axis.ticks.length = unit(0.1, "cm"),
          legend.position="none",
          strip.text.x = element_text(size = 12, colour = "black"),
          strip.background = element_rect(fill="gray88"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "ghostwhite"),
          strip.placement = 'outside',
          plot.caption = element_text(color = "Black", size = 10))
      
      
      png(file = paste("../../output/1_SNPs_identification/", folders[[z]], "/mp",  asso[[a]], ".png", sep = ""),
          width = 1200, height = 900)
      
      print(manhatan_plot)
      
      dev.off()
    
    }
  
  # Creates the QQplot:
  
  q_q_plot <- gg_qqplot(man_table$P) +
    theme_bw(base_size = 24) +
    annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1 + 0.15 * 3,
             label = sprintf("lambda = %.2f", inflation(man_table$P)), size = 8) +
    theme(axis.title = element_text(size = 25),
          axis.text.x  =  element_text(size = 25, colour = "black"),
          axis.text.y  = element_text(size = 25, hjust = 1, colour = "black"),
          legend.title = element_text(size = 25),
          legend.text  = element_text(size = 20),
          legend.key = element_rect(fill = "white"),
          legend.key.size = unit(1.3, "cm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin=unit(c(1, 3, 1, 3), "cm")) 

  
  png(file = paste("../../output/1_SNPs_identification/", folders[[z]], "/qq_plot_",  asso[[a]], ".png", sep = ""),
      width = 1200, height = 900)
  
  print(q_q_plot)
  
  dev.off()
  
}
  
}



