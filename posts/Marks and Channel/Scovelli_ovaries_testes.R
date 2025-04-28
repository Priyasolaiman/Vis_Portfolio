#Load in the two Libraries 
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(EnhancedVolcano) 

#Read in the count file 
S_ovaries <- read.csv("C:/Users/farli/OneDrive/Documents/sygnathus scovelli.rnaseq data/Final_count_of_all_tissue_syngnathus_scovelli/final_count_ss_ovaries.csv", header = TRUE, row.names = 1, sep = ",")
S_testes_oldref <- read.csv("C:/Users/farli/OneDrive/Documents/sygnathus scovelli.rnaseq data/Final_count_of_all_tissue_syngnathus_scovelli/final_count_ss_testes_oldrefgen.csv", header = TRUE, row.names = 1, sep =",")

#stretching the data in longer format as ggplot can't process multi column data at once
S_tst <- S_testes_oldref %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Counts")
S_ov <- S_ovaries %>% 
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Counts")

#Distribution plot across samples to test batch effects
#Testis, Old Reference genome(syngnathus scovelli)
ggplot(S_tst, aes(x = log10(Counts+1), fill = Sample)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  theme_minimal() +
  labs(title = "Gene Expression Distribution Across S_Testis", 
       x = "Log10(Gene Counts +1)", y = "Frequency")

#Ovaries, Syngnathus scovelli
ggplot(S_ov, aes(x = log10(Counts + 1), fill = Sample))+
  geom_histogram(alpha = 0.8, position = "identity", bins = 50)+
  theme_minimal()+
  labs(title = "Gene Expression Distribution Across S_Ovaries",
      x = "Log10(Gene Counts +1)", y = "Frequency")

#Histogram to check overall gene distribution
#making a dataframe of sum of all row count
S_tst_hist <- data.frame(rowSums(S_testes_oldref)) 
S_ov_hist <- data.frame(rowSums(S_ovaries))

#changing the column names
names(S_tst_hist)[1] <- "c1"
names(S_ov_hist)[1] <- "c3"

x_limits <- c(1, 1e07)  # Adjust based on your data range
y_limits <- c(0, 1000)   # Adjust to match both plots
y_breaks <- seq(0, 1000, by = 250)  # Ensure consistent y-axis breaks
x_breaks <- c(10, 1000, 100000, 1e07)  # Adjust log-scale

# Example for first dataset
ggplot(S_tst_hist, aes(x = c1)) + 
  geom_histogram(fill = "blue4", color = "black", bins = 100) +
  scale_x_log10(limits = x_limits, breaks = x_breaks) + 
  scale_y_continuous(limits = y_limits, breaks = y_breaks) + 
  labs(title = "Distribution of Overall Gene Counts Across S.scoveli Testis", 
       x = "Log(total counts)", y = "Frequency") +
  theme_gray()

ggplot(S_ov_hist, aes(x = c3)) + 
  geom_histogram(fill = "blue4", color = "black", bins = 100) +
  scale_x_log10(limits = x_limits, breaks = x_breaks) + 
  scale_y_continuous(limits = y_limits, breaks = y_breaks) + 
  labs(title = "Distribution of Overall Gene Counts Across S.scoveli Ovaries", 
       x = "Log(total counts)", y = "Frequency") +
  theme_gray()


#wanted to combine all the file together but different rows means different gene number
S_OT <- cbind(S_ovaries, S_testes_oldref)
#Filter out anything that has a sum of less than 10 but be careful regarding tissue
#S_OT <-  S_OT[rowSums(S_OT) > 10, ]
#Subset the Counts data for each of the different conditions 
All <- S_OT[, c(1:12)]
SFO_vs_SPT_count_table <- S_OT[, c(1:5, 8:12)]
SFO_vs_SNPT_count_table <- S_OT[, c(1:5, 6:7)]

##test
#Create the conditions for each of them 
All_condition <- c(rep("SFO",5), rep("SNPT",2), rep("SPT",5))
SFO_vs_SNPT_condition <- c(rep("SFO", 5), rep("SNPT", 2))
SFO_vs_SPT_condition <-  c(rep("SFO", 5), rep("SPT", 5))
###########################

#test
coldata_ALL <- data.frame(row.names = colnames(All), All_condition)
coldata_SFO_vs_SNPT <- data.frame(row.names = colnames(SFO_vs_SNPT_count_table), SFO_vs_SNPT_condition)
coldata_SFO_vs_SPT <- data.frame(row.names = colnames(SFO_vs_SPT_count_table), SFO_vs_SPT_condition)
############################

dds_ALL <- DESeqDataSetFromMatrix(countData = All, 
                                  colData = coldata_ALL, 
                                  design = ~All_condition)

dds_SFO_vs_SNPT <-  DESeqDataSetFromMatrix(countData = SFO_vs_SNPT_count_table, 
                                           colData = coldata_SFO_vs_SNPT,
                                           design = ~SFO_vs_SNPT_condition)
dds_SFO_vs_SPT <-  DESeqDataSetFromMatrix(countData = SFO_vs_SPT_count_table,
                                          colData = coldata_SFO_vs_SPT,
                                          design = ~SFO_vs_SPT_condition)

################################
dds_ALL <- DESeq(dds_ALL)
dds_SFO_vs_SNPT <- DESeq(dds_SFO_vs_SNPT)
dds_SFO_vs_SPT <- DESeq(dds_SFO_vs_SPT)

###########################
# Calling results without any arguments will extract the 
# estimated log2 fold changes and p values for the last variable in the design formula
#res_FO_vs_MT <- results(dds_FO_vs_MT)
res_all <- results(dds_ALL)
res_SFO_vs_SNPT <- results(dds_SFO_vs_SNPT)
res_SFO_vs_SPT <- results(dds_SFO_vs_SPT)
res_all

#mcols is basically shows metadata column names
mcols(res_all, use.names = TRUE)
#basemean=average of normalize count,log2FoldChange
#is the effect size estimate.lfcSE, 
#the standard error estimate for the log2 fold change estimate
#p value indicates the probability that a fold change as strong as the observed one, or even stronger,
sum(res_SFO_vs_SNPT$padj < 0.05, na.rm = TRUE)
sum(res_SFO_vs_SNPT$padj < 0.05, na.rm = TRUE)
###########################
#removing na values
sigs_all <- na.omit(res_all)
sigs_SFO_vs_SNPT <- na.omit(res_SFO_vs_SNPT)
sigs_SFO_vs_SPT <- na.omit(res_SFO_vs_SPT)


#sigs_FO_vs_MT <- sigs_FO_vs_MT[sigs_FO_vs_MT$padj < 0.05,]
sigs_SFO_vs_SNPT <- sigs_SFO_vs_SNPT[sigs_SFO_vs_SNPT$padj < 0.05,]
sigs_SFO_vs_SPT <- sigs_SFO_vs_SPT[sigs_SFO_vs_SPT$padj < 0.05,]
sigs_all <- sigs_all[sigs_all$padj < 0.05,]

#To see first few rows with lowest log2fold value/strongest down-regulation
head(sigs_SFO_vs_SNPT[order(sigs_SFO_vs_SNPT$log2FoldChange), ])
#To see strongest up-regulated genes
tail(sigs_SFO_vs_SNPT[order(sigs_SFO_vs_SNPT$log2FoldChange), ])

#plotting
#MA plot, x axis = mean expression & y axix = log2fold, c(-5,5) means log2fold value -5 to +5
#plotMA(res_all, ylim = c(-5,5))
plotMA(res_SFO_vs_SNPT, ylim = c(-5,5), main = "S_Ovaries vs Non Pregnant S_Testis")
plotMA(res_SFO_vs_SPT, ylim = c(-5,5), main = "S_Ovaries vs Pregnant S_Testis")

#plotDispEsts visualizes DESeq2 ’s dispersion estimates
#plotDispEsts(dds_ALL, ylim = c(1e-6, 1e1))
plotDispEsts(dds_SFO_vs_SNPT, ylim = c(1e-6, 1e1), main = "S_Ovaries vs Non Pregnant S_Testis")
plotDispEsts(dds_SFO_vs_SPT, ylim = c(1e-6, 1e1), main = "S_Ovaries vs Pregnant S_Testis")

#Histogram of P value
#hist(res_all$pvalue, breaks = 20, col = "grey")
hist(res_SFO_vs_SNPT$pvalue, 
     breaks = 20, 
     col = "grey", 
     main = "P-value Distribution S_Ov vs S_NPT", 
     xlab = "P-value", 
     ylab = "Frequency")

hist(res_SFO_vs_SPT$pvalue, 
     breaks = 20, 
     col = "grey", 
     main = "P-value Distribution S_ov vs S_PT", 
     xlab = "P-value", 
     ylab = "Frequency")


#histogram with adjusted p value
hist(res_SFO_vs_SNPT$padj, 
     breaks = 20, 
     col = "grey", 
     main = "P-value Distribution S_Ov vs S_NPT", 
     xlab = "P-value", 
     ylab = "Frequency")

hist(res_SFO_vs_SPT$padj, 
     breaks = 20, 
     col = "grey", 
     main = "P-value Distribution S_ov vs S_PT", 
     xlab = "P-value", 
     ylab = "Frequency")

# create bins using the quantile function  
qs <- c( 0, quantile( res_all$baseMean[res_all$baseMean > 0], 0:7/7 ) ) 
# "cut" the genes into the bins 
bins <- cut( res_all$baseMean, qs )  
# rename the levels of the bins using the middle point  
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)])) 
# calculate the ratio of £p£ values less than .01 for each bin  
ratios <- tapply( res_all$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) ) 
# plot these ratios  
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

attr(res_all,"filterThreshold") 
## 40% ## 165  
plot(attr(res_all,"filterNumRej"),
     type="b", xlab="quantiles of 'baseMean'", 
     ylab="number of rejections")

res_all$ge <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )

#write.csv(as.data.frame(res_all), "Result_all.csv", row.names = TRUE)
write.csv(as.data.frame(res_SFO_vs_SNPT), "Result_SFO_vs_SNPT.csv", row.names = TRUE)
write.csv(as.data.frame(res_SFO_vs_SPT), "Result_SFO_vs_SPT.csv", row.names = TRUE)

#rlog transform for application not for differntial testing
rld <- rlog(dds_ALL)

#PCAplot, plotPCA which comes with DESeq2.
# Run PCA and store the ggplot object
pca_plot <- plotPCA(rld, intgroup = "All_condition")

# Customize colors using ggplot2
pca_plot + scale_color_manual(values = c("#550000", "#AA0000", "#FF0000", 
                                         "#005500", "#00AA00", "#00FF00", 
                                         "#000055", "#0000AA", "#0000FF", 
                                         "#550055", "#AA00AA", "#FF00FF"))+
  labs(title = "PCA Plot of S_ovaries and S_Testis", 
       x = "PC1:84% Variance", 
       y = "PC2:4% Variance", 
       color = "Sample Groups") +  # Add legend title
  theme_minimal()

#Heatmap with gene clustering
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2(assay(rld)[topVarGenes, ], 
          scale = "row", trace = "none", 
          dendrogram = "column", 
          col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
          cexRow = 0.6,  # Adjust row label size
          cexCol = 0.7,
          main = "S_Ov vs S_T")  # Adjust column label size)

#venn diagram
#install.packages("VennDiagram")
#install.packages("ggVennDiagram")
library(VennDiagram)
library(ggVennDiagram)


df_SFO_vs_SNPT <- as.data.frame(sigs_SFO_vs_SNPT)
df_SFO_vs_SPT <- as.data.frame(sigs_SFO_vs_SPT)
df_all <-  as.data.frame(sigs_all)

###########################
write.csv(df_SFO_vs_SNPT, "df_SFO_SNPT.csv",row.names = TRUE)
write.csv(df_SFO_vs_SPT, "df_SFO_vs_SPT.csv", row.names = TRUE)

df.top_SFO_vs_SNPT <- df_SFO_vs_SNPT[df_SFO_vs_SNPT$baseMean > 50 & (abs(df_SFO_vs_SNPT$log2FoldChange) > 2.5),]
df.top_SFO_vs_SPT <- df_SFO_vs_SPT[df_SFO_vs_SPT$baseMean > 50 & (abs(df_SFO_vs_SPT$log2FoldChange) > 2.5),]
#df.top_L_vs_M <- df.top_L_vs_M[order(df.top_L_vs_M$log2FoldChange, decreasing = TRUE),]

#df.top_FO_vs_MT <- df.top_FO_vs_MT[order(df.top_FO_vs_MT$padj,df.top_FO_vs_MT$log2FoldChange,decreasing = c(TRUE, TRUE)),]
df.top_SFO_vs_SNPT <- df.top_SFO_vs_SNPT[order(df.top_SFO_vs_SNPT$padj, 
                                               df.top_SFO_vs_SNPT$log2FoldChange,
                                               decreasing = c(TRUE, TRUE)),]
                
df.top_SFO_vs_SPT <- df.top_SFO_vs_SPT[order(df.top_SFO_vs_SPT$padj, 
                                             df.top_SFO_vs_SPT$log2FoldChange,
                                             decreasing = c(TRUE, TRUE)),]

write.csv(df.top_SFO_vs_SNPT, "Scovelli_Female_Ovaries_vs_Nonpregnant_Male_Testis_DGE(0.05).csv", row.names = TRUE)
write.csv(df.top_SFO_vs_SPT, "Scovelli_Female_Ovaries_vs_Pregnant_Male_Testis_DGE(0.05)", row.names = TRUE)
#write.csv(df.top_L_vs_M, "Liver_vs_Muscle_DGE(0.01).csv", row.names = TRUE)

selected_gene_names_SFO_vs_SNPT <- row.names(head(df.top_SFO_vs_SNPT, 10))
selected_gene_names_SFO_vs_SPT <- row.names(head(df.top_SFO_vs_SPT, 10))
gene_names_SFO_vs_SNPT <- row.names(df_SFO_vs_SNPT)
gene_names_SFO_vs_SPT <- row.names(df_SFO_vs_SPT)


library(EnhancedVolcano) 

# Remove NA values in padj to avoid plotting issues
df_SFO_vs_SPT <- df_SFO_vs_SPT[!is.na(df_SFO_vs_SPT$padj), ]

# Assign colors based on log2FoldChange thresholds
OT_volcano <- ifelse(df_SFO_vs_SPT$log2FoldChange > 2.5, "orange3", #up in female ovaries
                     ifelse(df_SFO_vs_SPT$log2FoldChange < - 2.5, "maroon4", #down in female ovaries
                            "black"))

OT_volcano[is.na(OT_volcano)] <- "black"
names(OT_volcano)[OT_volcano == "orange3"] <- "Up FeMale_Ov"
names(OT_volcano)[OT_volcano== "black"] <- "NS"
names(OT_volcano)[OT_volcano == "maroon4"] <- "Up PregMale_T"

# Generate EnhancedVolcano plot
EnhancedVolcano(df_SFO_vs_SPT, 
                x = "log2FoldChange", 
                y = "padj", 
                lab = NA,  # Label genes
                pCutoff = 0.05, 
                title = "Female Ovaries vs Male Pregnant Testis", 
                FCcutoff = 2.5, 
                colCustom = OT_volcano,  # Assign custom colors
                colAlpha = 0.8,  # Adjust transparency
                drawConnectors = TRUE, 
                widthConnectors = 1.0, 
                boxedLabels = TRUE, 
                colConnectors = "black")

#ggsave("VP_L_vs_B.png", plot= VP_L_vs_B, dpi = 1000, width = 8, height = 8)
OT_volcano <- ifelse(df_SFO_vs_SPT$log2FoldChange < -2.5, "orange3",
                    ifelse(df_SFO_vs_SPT$log2FoldChange < 2.5, "maroon4",
                           "black"))
OT_volcano[is.na(OT_volcano)] <- "black"
names(OT_volcano)[OT_volcano == "orange3"] <- "Up Male_PT"
names(OT_volcano)[OT_volcano== "black"] <- "NS"
names(OT_volcano)[OT_volcano == "maroon4"] <- "Up Female_Ov"

EnhancedVolcano(df_SFO_vs_SPT, x = "log2FoldChange", y = "padj", lab = NA,
                pCutoff = 0.05, title = "Female Ovaries vs Male Pregnant Testis",
                FCcutoff = 2.5, colCustom = OT_volcano, colAlpha = 4/5,
                drawConnectors = TRUE, widthConnectors = 1.0, boxedLabels = TRUE,
                colConnectors = "black")


EnhancedVolcano(df_FO_vs_MT, x = "log2FoldChange", y = "padj", lab = gene_names_FO_vs_MT , 
                selectLab = c('LOC133154020', 'dmrt1'),
                pCutoff = 0.05, title = "Ovary vs Testes", 
                FCcutoff = 2.5, 
                col = c('black', 'black', 'springgreen3', 'aquamarine4')
                ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                , boxedLabels = TRUE, colConnectors = 'black')

keyvals2 <- ifelse(
  df_B_vs_M$log2FoldChange < -2.5, 'coral3',  
  ifelse(df_B_vs_M$log2FoldChange > 2.5, 'darkred',  
         'black'))
keyvals2[is.na(keyvals2)] <- 'black'
names(keyvals2)[keyvals2 == 'darkred'] <- 'Up in Muscle'
names(keyvals2)[keyvals2 == 'black'] <- 'NS'
names(keyvals2)[keyvals2 == 'coral3'] <- 'Down in Muscle'

B_v_M <- EnhancedVolcano(df_B_vs_M, x = "log2FoldChange", y = "padj", lab = NA,
                         pCutoff = 0.05, title = "B) Brain vs Muscle", 
                         FCcutoff = 2.5, 
                         subtitle = NULL,
                         colCustom = keyvals2
                         ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                         , boxedLabels = TRUE, colConnectors = 'black')# +coord_flip()

keyvals3 <- ifelse(
  df_L_vs_M$log2FoldChange < -2.5, 'navyblue',  
  ifelse(df_L_vs_M$log2FoldChange > 2.5, 'darkred',  
         'black'))
keyvals2[is.na(keyvals3)] <- 'black'
names(keyvals3)[keyvals3 == 'darkred'] <- 'Up in Liver'
names(keyvals3)[keyvals3 == 'black'] <- 'NS'
names(keyvals3)[keyvals3 == 'navyblue'] <- 'Down in Liver'

L_v_M <- EnhancedVolcano(df_L_vs_M, x = "log2FoldChange", y = "padj", lab = NA , 
                         pCutoff = 0.05, title = "C) Liver vs Muscle",
                         subtitle = NULL,
                         FCcutoff = 2, 
                         colCustom = keyvals3
                         ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                         , boxedLabels = TRUE, colConnectors = 'black')# + coord_flip()

library(ggpubr)
ggarrange(B_V_L, B_v_M, L_v_M, ncol = 3)


EnhancedVolcano(df_B_vs_LM, x = "log2FoldChange", y = "padj", lab = gene_names_B_vs_LM , 
                pCutoff = 0.05, title = "Brain vs Liver & Muscle",
                FCcutoff = 2, 
                col = c('black', 'black', 'springgreen3', 'aquamarine4')
                ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                , boxedLabels = TRUE, colConnectors = 'black')

EnhancedVolcano(df_L_vs_BM, x = "log2FoldChange", y = "padj", lab = gene_names_L_vs_BM , 
                pCutoff = 0.05, title = "Liver vs Brain & Muscle",
                FCcutoff = 2, 
                col = c('black', 'black', 'springgreen3', 'aquamarine4')
                ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                , boxedLabels = TRUE, colConnectors = 'black')

EnhancedVolcano(df_M_vs_BL, x = "log2FoldChange", y = "padj", lab = gene_names_M_vs_BL , 
                pCutoff = 0.05, title = "Muscle vs Brain & Liver",
                FCcutoff = 2, 
                col = c('black', 'black', 'springgreen3', 'aquamarine4')
                ,colAlpha = 4/5, drawConnectors = TRUE, widthConnectors = 1.0
                , boxedLabels = TRUE, colConnectors = 'black')





#normalized_counts_B_v_L <- counts(dds_B_vs_L, normalized=TRUE)

