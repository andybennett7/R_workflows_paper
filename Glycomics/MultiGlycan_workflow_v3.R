############# MultiGlycan Workflow #######

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

## May need to set wd to where files are
setwd("C:\\Users\\Andy Bennett\\Desktop\\R\\R code manuscript\\Example glycan data")
fileNames <- Sys.glob("*.csv")
samples <- lapply(fileNames, function(x) {
  read.csv(x)})

# Keep only "All Adducts"
for (fileName in fileNames) {
  x <- read.csv(fileName)
  all_glycans <- separate(x, col=Glycan, into=c('Glycan', 'Adducts', 'X'), sep=' ', remove=TRUE)
  all_glycans <- separate(all_glycans, col=Glycan, into=c('HexNAc', 'Hex','Fuc', 'NeuAc', 'NeuGc'), sep='-', remove=FALSE)
  
  all_adducts <- all_glycans[!grepl("Ratio", all_glycans$X),]
  all_adducts <- all_adducts[!grepl("Protonated", all_adducts$Adducts),]
  all_adducts <- within(all_adducts, rm("X"))
  all_adducts <- within(all_adducts, rm("Adducts"))
  
  #### remove columns with no data
  # change N/A to NA
  is.na(all_adducts) <- all_adducts == "N/A"
  
  # Eliminate columns with no data
  all_adducts <-  all_adducts[,!sapply(all_adducts, function(x) mean(is.na(x)))>0.95]
  
  #name data frame
  new_name <- sub(".csv", "", fileName)
  #regroup into list "samples_2"
  write.csv(all_adducts, file = paste0("cleaned_",new_name,".csv"), row.names = FALSE)
}
#sample files
fileNames_2 <- Sys.glob("cleaned_*")
samples_2 <- lapply(fileNames_2, function(x) {
  read.csv(x)})

#Create Master_Glycan_list 
Glycan_master_list <- lapply(fileNames_2, function(x) {
  read.csv(x)["Glycan"]})

# unlist to combine into a vector
Glycan_master_list <- data.frame(unlist(Glycan_master_list))
Glycan_master_list_cleaned <- unique(Glycan_master_list)
rownames(Glycan_master_list_cleaned) = NULL
colnames(Glycan_master_list_cleaned) <- "Glycan"

# Create a new data frame to store the aligned data
aligned_df <- Glycan_master_list_cleaned

# Align data by GP master list, now renamed as aligned_df
for (df in samples_2) {
  aligned_df <- merge(aligned_df, df, all = TRUE)
}

# remove rows with NAs after alignment 8 was used because that was the number of
# samples. should be adjusted for other data sets
aligned_df <- aligned_df[rowSums(is.na(aligned_df)) < 12, ]

# reorder sample columns
aligned_df2 <- aligned_df[,7:ncol(aligned_df)] %>% 
  select(sort(names(.)))
# bind to glycan
aligned_df2 <- cbind(aligned_df[1:6], aligned_df2)

#normalize based on abundance
smpl_start <- 7
column_sums <- colSums(aligned_df2[smpl_start:ncol(aligned_df2)], na.rm = TRUE)

# Divide each analyte by sum of that sample 
a <- function(x){x/column_sums}
norm_aligned <- apply(aligned_df2[smpl_start:ncol(aligned_df2)],2,a)
norm_aligned <- as.data.frame(norm_aligned)
#add back the columns lost in calculations
info_column_end <- (as.numeric(unlist(smpl_start)) - 1)
norm_aligned <- cbind(aligned_df2[1:info_column_end], norm_aligned)


# remove all rows with any NAs. Important to remove AFTER normalization for more
# accurate results
norm_aligned <- na.omit(norm_aligned)

####### STATS #######
#Fold change

data_AD <- norm_aligned %>% select(contains(c("AD")))
data_C <- norm_aligned %>% select(matches(c("C[1-9]")))
norm_aligned$AD_avg <- rowMeans(data_AD)
norm_aligned$Control_avg <- rowMeans(data_C)

norm_aligned$AD_vs_C_fold_change <- (norm_aligned$AD_avg)/(norm_aligned$Control_avg)

#Log2 calculation
norm_aligned$log2fc <- log(norm_aligned$AD_vs_C_fold_change, 2)


#pval Welch's t.test
# AD vs C
norm_aligned$AD_C_pval <- sapply(1:nrow(norm_aligned), function(i) t.test(as.numeric(as.character(unlist(data_AD[i,1:6]))), as.numeric(as.character(unlist(data_C[i,1:6]))))[c("p.value")])
norm_aligned$AD_C_pval <- as.numeric(norm_aligned$AD_C_pval)
#sig_before <- sum(norm_aligned$AD_C_pval < 0.05)

####### SUBSET BY TYPE #######

High_Mannose <- subset(norm_aligned, norm_aligned$HexNAc < 3 & norm_aligned$Fuc==0 & norm_aligned$NeuAc==0 & norm_aligned$NeuGc==0)
Sialylated_NeuAc <- subset(norm_aligned, norm_aligned$HexNAc > 2 & norm_aligned$Fuc==0 & norm_aligned$NeuAc>0 & norm_aligned$NeuGc==0)
Sialylated_NeuGc <- subset(norm_aligned, norm_aligned$HexNAc > 2 & norm_aligned$Fuc==0 & norm_aligned$NeuAc==0 & norm_aligned$NeuGc>0)
Fucosylated <- subset(norm_aligned, norm_aligned$HexNAc > 1 & norm_aligned$Fuc>0 & norm_aligned$NeuAc==0 & norm_aligned$NeuGc==0)
SialoFucosylated <-  subset(norm_aligned, norm_aligned$HexNAc > 2 & norm_aligned$Fuc>=1 & (norm_aligned$NeuAc>=1 | norm_aligned$NeuGc>=1))
Others <- subset(norm_aligned, norm_aligned$HexNAc > 2 & norm_aligned$Fuc==0 & norm_aligned$NeuAc==0 & norm_aligned$NeuGc==0)

#check if subset groups is equal to total glycans
a <- c(High_Mannose$Glycan,Sialylated_NeuAc$Glycan,Sialylated_NeuGc$Glycan,SialoFucosylated $Glycan,Fucosylated$Glycan,Others$Glycan)
b <- norm_aligned$Glycan

#Should be TRUE! #check that glycans have been sorted correctly
length(a)==length(b)


####### VISUALIZATION #######
# Glycan colors (SNFG from `Essentials of Glycobiology`)
color_code <- data.frame(matrix(NA, nrow = 10, ncol = 4))

colnames(color_code) <- c("Color", "CMYK", "RBG", "HEX")
color_code$Color <- c("White","Blue","Green","Yellow","Light_blue","Pink","Magenta","Brown","Orange","Red")
color_code$CMYK <- c("0/0/0/0", "100/50/0/0", "100/0/100/0", "0/15/100/0", "41/5/3/0", "0/47/24/0", "38/88/0/0", "32/48/76/13", "0/65/100/0", "0/100/100/0")
color_code$RBG <- c("255/255/255", "0/114/188","0/166/81", "255/212/0", "143/204/233", "246/158/161", "165/67/153", "161/122/77", "244/121/32", "237/28/36")

color_code$HEX <- c("#FFFFFF", "#0072BC", "#00A651", "#FFD400", "#8FCCE9", "#F69EA1", "#A54399", "#A17A4D", "#F47920", "#ED1C24")

# 2D PCA plot
#packages

norm_pca <- norm_aligned[,c(7:18)] 

norm_pca <- data.frame(t(norm_pca))
norm_pca$Sample <- c("AD (N=6)", "AD (N=6)","AD (N=6)","AD (N=6)","AD (N=6)","AD (N=6)",  "C (N=6)", "C (N=6)","C (N=6)","C (N=6)","C (N=6)","C (N=6)" )
norm_pca <- norm_pca %>% relocate(Sample)
data.pca <- prcomp(norm_pca[,c(2:ncol(norm_pca))], 
                   center = TRUE, 
                   scale. = TRUE)

fviz_pca_ind(data.pca,
             geom.ind = c("point"), # show points only (nbut not "text")
             col.ind = norm_pca$Sample, # color by groups
             palette = c("purple2", "orange2"),
             legend.title = "Groups", 
             title = "PCA of Glycans in AD vs C", 
             mean.point = F)

ggsave(here("Example glycan data","plots", "pca_glycans.svg"))

#Scree plot
fviz_eig(data.pca, addlabels = TRUE)
ggsave(here("Example glycan data","plots", "scree_glycans.svg"))



# Volcano

#select data
clean_data <- norm_aligned

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

#Total Glycans
total_gly <- (ns+down_reg+up_reg)

ymax <- -log(min(unlist(clean_data$AD_C_pval)), 10)
xmin <- min(clean_data$log2fc)
xmax <- max(clean_data$log2fc)

annotation <- data.frame(
  x = c((xmin*0.75),(xmax*0.75), -Inf),
  y = c((ymax*0.5), (ymax*0.5), Inf),
  label = c(down_reg, up_reg, paste0("Total Glycans: ", total_gly)))

#Theme Set
theme_set(theme_bw())

v <- ggplot(data = clean_data, aes(x = log2fc, y = -log10(unlist(AD_C_pval)), col = diffexpressed, )) + 
  #label = data_label
  
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(1.3), col = "black", linetype = "dashed", linewidth = 1) +
  geom_point(size = 1.5, alpha = 0.7) + 
  #geom_text_repel(max.overlaps = Inf) +
  
  geom_label(data=annotation, aes(x=x, y=y, label=label, hjust="inward", vjust="inward"),
             color="black", 
             size=4, angle=0, fontface="bold") + 
  ggtitle(graph_title) + 
  labs(col = "Glycans") + labs(x = expression(log[2]*"fold change"), y = expression(-log[10]*"p-value")) +
  
  scale_color_manual(values = c("firebrick", "grey60", "forestgreen"), 
                     labels = c("Downregulated", "Outside Statistical Parameters", "Upregulated"))

### optional colors for volcano:
### down(purple3, red3, blue4, firebrick)
### up (yellow3, green3, orange3, forestgreen)

#p+theme(axis.text = element_text(size = 14))
v <- v+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 12))+theme(legend.text = element_text(size = 10))+guides(color = guide_legend(override.aes = list(size = 4)))+theme(legend.title = element_blank())+theme(legend.position = "bottom")+theme(plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

v

ggsave(here("Example glycan data","plots", "volcano_glycans.svg"))

#heatmap data set up
heat_data <- subset(norm_aligned, norm_aligned$AD_C_pval < 0.05)
heat_data <- heat_data[,c(1,7:18)] 
#remove row names
glycans <- heat_data$Glycan
glycans <- data.frame(glycans)
heat_data2 <- heat_data[,-1]
heat_map_data <- as.matrix(heat_data2)
row.names(heat_map_data) <- heat_data$Glycan

colr <- colorRampPalette(c("dodgerblue3", "gray0", "green3"))(100)


# heatmap
h <- pheatmap(
  heat_map_data  ,
  color = colr,
  scale = "row",
  angle_col = "0",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "AD vs C",
  fontsize = 10, 
  show_rownames = TRUE,
  fontsize_row = 8
)

ggsave(here("Example glycan data","plots", "heatmap_glycans.svg"), plot = h )

#composition by glycan type
composition <- data.frame(matrix(NA, nrow = 5, ncol = 2))

colnames(composition) <- c("AD", "C")

composition$Glycan_Type <- c("High Mannose", "Fucosylated", "Sialylated_NeuAc", "SialoFucosylated", "Others" )
# Note: NeuGc was not included, but can easily be added for other experiments/data
composition$AD <- c(sum(High_Mannose$AD_avg), sum(Fucosylated$AD_avg), sum(Sialylated_NeuAc$AD_avg), sum(SialoFucosylated$AD_avg), sum(Others$AD_avg))
composition$C <- c(sum(High_Mannose$Control_avg), sum(Fucosylated$Control_avg), sum(Sialylated_NeuAc$Control_avg), sum(SialoFucosylated$Control_avg), sum(Others$Control_avg))

composition$AD <- composition$AD*100
composition$C <- composition$C*100
AD_sum_glycans <- sum(composition$AD)
C_sum_glycans <- sum(composition$C)
composition$AD <- composition$AD/AD_sum_glycans
composition$C <- composition$C/C_sum_glycans

glycan_colors <- c("#00A651","#ED1C24","#A54399", "#FFD400",'gray30')

comp_3 <- pivot_longer(composition, cols = c("AD", "C"), names_to = "Group", values_to = "Percentage")
comp_3$Percentage <- round(comp_3$Percentage, 3)

ggplot(comp_3, aes(y = Percentage, x = Glycan_Type, fill = Group)) +
  geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values=c("purple2","orange2")) + ggtitle("Glycan Composition by Disease State") +theme(plot.title=element_text(hjust=0.5)) + geom_text(aes(label=Percentage), position=position_dodge(width=0.9), vjust=-0.25)

ggsave(here("Example glycan data","plots", "comp_by_disease.svg"))

ggplot(comp_3, aes(y = Percentage, x = Group, fill = Glycan_Type)) +
  geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values=c("#ED1C24","#00A651","gray30","#FFD400", "#A54399")) + ggtitle("Glycan Composition by Glycan Type") +theme(plot.title=element_text(hjust=0.5)) + geom_text(aes(label=Percentage), position=position_dodge(width=0.9), vjust=-0.25) 

ggsave(here("Example glycan data","plots", "comp_by__glycans.svg"))


#box plots (top 10 sig)
top_10_sig <- subset(norm_aligned, ((norm_aligned$log2fc < -1 | norm_aligned$log2fc > 1)) & norm_aligned$AD_C_pval < 0.05)
top_10_sig2 <- head(arrange(top_10_sig, AD_C_pval), 10)

top_10_list <- top_10_sig2$Glycan

for (glycan in top_10_list) { 
  column_number <- which(top_10_sig2$Glycan == glycan)
  data_subset <- top_10_sig2[column_number, 7:18]
  
  long_data <- data_subset %>%
    gather(key = "Groups", value = "value", 1:12) %>%
    mutate(group = case_when(
      row_number() <= 6 ~ "AD",
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
      scale_fill_manual(values = c("AD" = "purple2", "Control" = "orange2")) + geom_signif(comparisons = list(c("AD", "Control")), map_signif_level = TRUE)
    
    #stat_compare_means(comparisons = list(c("AD", "Control")), label = "p.signif",label.x = 1.5) +
    
    theme(axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 13),
          panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # Increase plot border linewidth here
    ) +
      ggtitle(paste0(glycan)) +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #axis.line = element_line(color = "black", size = 1.1),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")# Adds solid black lines for x and y axes
      )# This will remove the panel background as well
  }
  
  ggsave(here("Example glycan data","plots", paste0("boxplot_",glycan, ".svg")))
}  

# Export Results to Excel
write.xlsx(aligned_df2, "R_glycan_results.xlsx", sheetName = "original_data")
write.xlsx(norm_aligned, "R_glycan_results.xlsx", sheetName="stats", append=TRUE)
write.xlsx(top_10_sig, "R_glycan_results.xlsx", sheetName="significant", append=TRUE)
write.xlsx(High_Mannose, "R_glycan_results.xlsx", sheetName="High Mannose", append=TRUE)
write.xlsx(Fucosylated, "R_glycan_results.xlsx", sheetName="Fucosylated", append=TRUE)
write.xlsx(Sialylated_NeuAc, "R_glycan_results.xlsx", sheetName="Sialylated_NeuAc", append=TRUE)
write.xlsx(Sialylated_NeuGc, "R_glycan_results.xlsx", sheetName="Sialylated_NeuGc", append=TRUE)
write.xlsx(SialoFucosylated, "R_glycan_results.xlsx", sheetName="SialoFucosylated", append=TRUE)
write.xlsx(Others, "R_glycan_results.xlsx", sheetName="Others", append=TRUE)

