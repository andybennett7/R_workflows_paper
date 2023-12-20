##### Updated script to process  files from LipidSearch software.

#Load Library 
{
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(xlsx)
  library(pheatmap)
  library(ggpubr)
  library(ggsignif)
  library(here)
  library(ggfortify)
  library(FactoMineR)
  library(factoextra)
  library(corrr)
  library(ggcorrplot)
}

setwd(here("Example_lipid_data/"))
fileNames <- Sys.glob("*.txt")
for (fileName in fileNames) {
  #read in file
  data <- read.table(fileName, header=TRUE, skip=5, sep="\t")
  
  #removes low level and less confident results from software results
  data2 <- subset(data, data$Rt>2 & data$Area>1000 & data$t.Score<0.5)
  rownames(data2) <- NULL
  data2 <- arrange(data2, Rt)
  
  #keep only selected columns. can be adjusted according to project
  data3 <- data2 %>% select(-contains(c("Formula", "It.", "Pol.", "Scan", "CalcMz", "Delta", "TopRT", "PeakQuality", "QMethod", "Height", "Hwhm", "QuantInfo")))
  
  new_name <- sub("*.txt", "", fileName)
  
  #chose to write to Excel so multiople sheets can be used, as opposed to CSV file
  write.xlsx(data3, paste0("cleaned_",new_name,".xlsx"), sheetName = "filtered", row.names = FALSE)
  
  data4 <- data3 %>% 
    group_by(Class) %>% 
    summarise(Area = sum(Area))
  colnames(data4) <- c("Class",new_name)
  rownames(data4) <- NULL
  write.xlsx(data4 ,paste0("cleaned_",new_name,".xlsx"), sheetName = "Area", append = TRUE)
}

#Master dataframe with Class
aligned_df <- data.frame(c("CL","DLCL","MLCL","LPA","PA","LPC","PC","LPE","PE","LPG","PG","LPI","PI","LPS","PS","PIP","PIP2","PIP3","Cer","CerP","CerPE","Hex1SPH","Hex1Cer","Hex2Cer","Hex3Cer","CerG2GNAc1","CerG3GNAc1","CerG3GNAc2","GM3","GM2","GM1","GD1a","GD1b","GD2","GD3","GT1a","GT1b","GT1c","GT2","GT3","GQ1c","GQ1b","LSM","phSM","SM","SPH","SPHP","ST","AcHexSiE","AcHexStE","AcHexZyE","AcHexCmE","AcHexChE","ChE","D7ChE","CmE","DG","D5DG","MG","SiE","StE","TG","D5TG","ZyE","AcCa","AEA","Co","cPA","FA","LPEt","LPMe","OAHFA","PAF","PEt","PMe","SL","WE"))
colnames(aligned_df) <- "Class"
aligned_df <- data.frame(sort(aligned_df$Class))
colnames(aligned_df) <- "Class"
master_list <- aligned_df


fileNames_2 <- Sys.glob("cleaned*")
samples <- lapply(fileNames_2, function(x) {
  read_excel(x, sheet = "Area")})

# Align data by Class master list, now renamed as aligned_sum
for (df in samples) {
  #df <- df[,-1]
  aligned_df <- merge(aligned_df, df, by = "Class", all = TRUE)
}

aligned_df <- aligned_df %>% select(-contains(c("1.x", "1.y", "...1")))

#Excel file with aligned results
write.xlsx(aligned_df,"Example_lipid_data_results.xlsx", sheetName = "processed", row.names = FALSE)

#normalize based on abundance
smpl_start <- readline(prompt = "Column of first sample?: " )
column_sums <- colSums(aligned_df[smpl_start:ncol(aligned_df)], na.rm = TRUE)


# Divide each analyte by sum of that sample 
a <- function(x){x/column_sums}
norm_aligned <- apply(aligned_df[smpl_start:ncol(aligned_df)],1,a)
norm_aligned <- as.data.frame(norm_aligned)
norm_aligned <- t(norm_aligned)

#add back the columns lost in calculations
info_column_end <- (as.numeric(unlist(smpl_start)) - 1)
norm_aligned <- cbind(aligned_df[1:info_column_end], norm_aligned)

# reorder sample columns
norm_aligned2 <- norm_aligned[,2:ncol(norm_aligned)] %>% 
  select(sort(names(.)))
# bind to Lipid list
norm_aligned2 <- cbind(norm_aligned[1], norm_aligned2)
row.names(norm_aligned2) <-  NULL

# remove all rows with any NAs. Important to remove AFTER normalization for more accurate results
norm_aligned2 <- na.omit(norm_aligned2)

norm_2 <- pivot_longer(norm_aligned2, cols = -c(Class), names_to = "Sample", values_to = "Relative Abundance")


####### STATS #######
#Fold change

data_AD <- norm_aligned2 %>% select(contains(c("AD")))
data_C <- norm_aligned2 %>% select(matches(c("C[1-9]")))
norm_aligned2$AD_avg <- rowMeans(data_AD)
norm_aligned2$Control_avg <- rowMeans(data_C)

norm_aligned2$AD_vs_C_fold_change <- (norm_aligned2$AD_avg)/(norm_aligned2$Control_avg)

#Log2 calculation
norm_aligned2$log2fc <- log(norm_aligned2$AD_vs_C_fold_change, 2)


#pval Welch's t.test
# AD vs C
norm_aligned2$AD_C_pval <- sapply(1:nrow(norm_aligned2), function(i) t.test(as.numeric(as.character(unlist(data_AD[i,1:5]))), as.numeric(as.character(unlist(data_C[i,1:5]))))[c("p.value")])
norm_aligned2$AD_C_pval <- as.numeric(norm_aligned2$AD_C_pval)
#sig_before <- sum(norm_aligned$AD_C_pval < 0.05)

##### Visualization #####
# 2D PCA plot
#packages

norm_pca <- norm_aligned2[,c(2:11)] 

norm_pca <- data.frame(t(norm_pca))
norm_pca$Sample <- c("AD (N=5)", "AD (N=5)","AD (N=5)","AD (N=5)","AD (N=5)",  "C (N=5)", "C (N=5)","C (N=5)","C (N=5)","C (N=5)")
norm_pca <- norm_pca %>% relocate(Sample)
data.pca <- prcomp(norm_pca[,c(2:ncol(norm_pca))], 
                   center = TRUE, 
                   scale. = TRUE)

fviz_pca_ind(data.pca,
             geom.ind = c("point"), # show points only (nbut not "text")
             col.ind = norm_pca$Sample, # color by groups
             palette = c("purple2", "orange2"),
             legend.title = "Groups", 
             addEllipses = T,
             title = "PCA of Lipids in AD vs C", 
             mean.point = F
             )

ggsave(here("Example_lipid_data","plots", "pca_lipids.svg"))

#Scree plot
fviz_eig(data.pca, addlabels = TRUE)
ggsave(here("Example_lipid_data","plots", "scree_lipids.svg"))
?fviz_eig

# Volcano

#select data
clean_data <- norm_aligned2

#User Input Customization
graph_title <-  as.character(readline(prompt = "Graph Title?: " ))

{
  clean_data$diffexpressed <- "NO"
  clean_data$diffexpressed[clean_data$log2fc > 1 & clean_data$AD_C_pval < 0.05] <- "UP"
  clean_data$diffexpressed[clean_data$log2fc < -1 & clean_data$AD_C_pval < 0.05] <- "DOWN"
  
  up_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "UP")))
  down_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "DOWN")))
  ns <- as.numeric(sum(str_count(clean_data$diffexpressed, "NO")))
}

#Total Lipids
total_lipids <- (ns+down_reg+up_reg)

ymax <- -log(min(unlist(clean_data$AD_C_pval)), 10)
xmin <- min(clean_data$log2fc)
xmax <- max(clean_data$log2fc)

annotation <- data.frame(
  x = c((xmin*0.75),(xmax*0.75), -Inf),
  y = c((ymax*0.5), (ymax*0.5), Inf),
  label = c(down_reg, up_reg, paste0("Lipid Classes: ", total_lipids)))

#Theme Set
theme_set(theme_bw())

v <- ggplot(data = clean_data, aes(x = log2fc, y = -log10(unlist(AD_C_pval)), col = diffexpressed, )) + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(1.3), col = "black", linetype = "dashed", linewidth = 1) +
  geom_point(size = 1.5) + 
  geom_text_repel(label = norm_aligned2$Class, max.overlaps = Inf) +
  
  geom_label(data=annotation, aes(x=x, y=y, label=label, hjust="inward", vjust="inward"),
             color="black", 
             size=4, angle=0, fontface="bold") + 
  ggtitle(graph_title) + 
  labs(col = "Glycans") + labs(x = expression(log[2]*"fold change"), y = expression(-log[10]*"p-value")) +
  
  scale_color_manual(values = c( "grey60", "forestgreen"), 
                     labels = c( "Outside Statistical Parameters", "Upregulated"))

### optional colors for volcano:
### down(purple3, red3, blue4, firebrick)
### up (yellow3, green3, orange3, forestgreen)

v <- v+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 12))+theme(legend.text = element_text(size = 10))+guides(color = guide_legend(override.aes = list(size = 4)))+theme(legend.title = element_blank())+theme(legend.position = "bottom")+theme(plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

v

ggsave(here("Example_lipid_data","plots", "volcano_lipids.svg"))

#heatmap data set up
heat_data <- subset(norm_aligned2, norm_aligned2$AD_C_pval < 0.05)
heat_data <- heat_data[,c(1:11)] 
#remove row names
lipids <- heat_data$Class
lipids <- data.frame(lipids)
heat_data2 <- heat_data[,-1]
heat_map_data <- as.matrix(heat_data2)
row.names(heat_map_data) <- heat_data$Class

colr <- colorRampPalette(c("dodgerblue3", "gray0", "green3"))(100)


# heatmap
h <- pheatmap(
  heat_map_data  ,
  color = colr,
  scale = "row",
  angle_col = "0",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Heatmap of Statistically Significant Lipids in AD vs C",
  fontsize = 10, 
  show_rownames = TRUE,
  fontsize_row = 8
)
ggsave(here("Example_lipid_data","plots", "heatmap_lipids.svg"), plot = h)

#composition by lipid type
composition <- data.frame(matrix(NA, nrow = nrow(norm_aligned2), ncol = 2))

colnames(composition) <- c("AD", "C")

composition$Lipid_Type <- c(norm_aligned2$Class)
# Note: NeuGc was not included, but can easily be added for other experiments/data
composition$AD <- c(norm_aligned2$AD_avg)
composition$C <- c(norm_aligned2$Control_avg)

composition$AD <- composition$AD*100
composition$C <- composition$C*100
AD_sum_lipids <- sum(composition$AD)
C_sum_lipids <- sum(composition$C)
composition$AD <- composition$AD/AD_sum_lipids
composition$C <- composition$C/C_sum_lipids

comp_3 <- pivot_longer(composition, cols = c("AD", "C"), names_to = "Group", values_to = "Percentage")
comp_3$Percentage <- as.numeric(round(comp_3$Percentage, 3))

comp_3 <- arrange(comp_3, "Percentage")

comp_4 <- slice_max(comp_3, order_by = Percentage, n = 8)

comp_AD <- subset(comp_3, comp_3$Group == "AD")
comp_AD2 <-  slice_max(comp_AD, order_by = Percentage, n = 5) 
top5_AD <- comp_AD2$Lipid_Type

comp_C <- subset(comp_3, comp_3$Group == "C")
comp_C2 <-  slice_max(comp_C, order_by = Percentage, n = 5) 
top5_C <- comp_C2$Lipid_Type

comp_bar_AD <- filter(norm_2, Class %in% top5_AD)
comp_bar_AD <-  subset(comp_bar_AD, comp_bar_AD$Type == "AD")

comp_bar_C <- filter(norm_2, Class %in% top5_C)
comp_bar_C <-  subset(comp_bar_C, comp_bar_C$Type == "C")

bar_data <- rbind(comp_bar_AD, comp_bar_C)

theme_set(theme_bw())
ggplot(comp_4, aes(y = Percentage, x = Lipid_Type, fill = Group)) +
  geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values=c("purple2","orange2")) + ggtitle("Lipid Composition by Disease State (Top 4 Most Abundant)") +theme(plot.title=element_text(hjust=0.5)) + geom_text(aes(label=Percentage), position=position_dodge(width=0.9), vjust=-0.25)

ggsave(here("Example_lipid_data","plots", "comp_by_disease.svg"))

ggplot(bar_data, aes(y = `Relative Abundance`, x = Type, fill = Class)) +
  stat_summary(fun="mean",geom="col",position="dodge") + stat_summary(fun.data="mean_se",geom="errorbar", position = "dodge") + ggtitle("Top 5 Lipids by Lipid Type") +theme(plot.title=element_text(hjust=0.5))  

ggsave(here("Example_lipid_data","plots", "comp_by_lipids_error_bars.svg"))


#box plots (top 10 sig)
top_10_sig <- subset(norm_aligned2, ((norm_aligned2$log2fc < -1 | norm_aligned2$log2fc > 1)) & norm_aligned2$AD_C_pval < 0.05)
top_10_sig2 <- head(arrange(top_10_sig, AD_C_pval), 10)

top_10_list <- top_10_sig2$Class

for (lipid in top_10_list) { 
  column_number <- which(top_10_sig2$Class == lipid)
  data_subset <- top_10_sig2[column_number, 2:11]
  
  long_data <- data_subset %>%
    gather(key = "Groups", value = "value", 1:10) %>%
    mutate(group = case_when(
      row_number() <= 5 ~ "AD",
      TRUE ~ "Control"
    ))
  
  #max_value <- max(long_data$value, na.rm = TRUE)
  #y_limit <- max_value * 1.4
  theme_set(theme_bw())
  # Plotting the BoxPlot
  {
    ggplot(long_data, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha=0.7, width = 0.4, lwd = 0.4, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0.7,
                   outlier.color='black') +
      xlab("Groups") +
      ylab("Relative Abundance (%)") +
      scale_fill_manual(values = c("AD" = "purple2", "Control" = "orange2")) + #geom_signif(comparisons = list(c("AD", "Control")), map_signif_level = TRUE)
      
      stat_compare_means(comparisons = list(c("AD", "Control")), label = "p.signif",label.x = 1.5) +
      
      theme(axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 13),
            panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # Increase plot border linewidth here
      ) +
      ggtitle(lipid) +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #axis.line = element_line(color = "black", size = 1.1),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")# Adds solid black lines for x and y axes
      )# This will remove the panel background as well
  }
  ?stat_compare_means
  ggsave(here("Example_lipid_data","plots", paste0("boxplot_",lipid, ".svg")))
}  

# Export Results to Excel
write.xlsx(aligned_df, "R_lipid_results.xlsx", sheetName = "grouped area")
write.xlsx(norm_aligned2, "R_lipid_results.xlsx", sheetName="stats", append=TRUE)
write.xlsx(top_10_sig, "R_lipid_results.xlsx", sheetName="significant", append=TRUE)
write.xlsx(clean_data, "R_lipid_results.xlsx", sheetName="found in all", append=TRUE)
write.xlsx(heat_data, "R_lipid_results.xlsx", sheetName="heatmap", append=TRUE)
write.xlsx(composition, "R_lipid_results.xlsx", sheetName="composition", append=TRUE)
