
library(tidyverse)
library(here)
library(readxl)
library(xlsx)

#setwd(here())
fileNames <- Sys.glob("*.xlsx")
fileNames
for (fileName in fileNames) {
  data <-  as.data.frame((read_xlsx(path = fileName, sheet = "Spectra", col_types = "text")))
  colnames(data) <- gsub("\n", "_", colnames(data))
  data2 <- subset(data, !is.na(data$Glycans_NHFAGNa))
  #data$Score > 300
  #combine columns into GP ID
  #GPs = protein + glycan + modification + z(charge + scan time
  data3 <- mutate(data2, GPs = paste(data2[, 3], data2[, 4],  data2[, 7], sep = "_"))
  #data2[,5], this is for modification
  new_name <- sub(".xlsx", "", fileName)
  data3$Sample <- paste(new_name)
  
  #name data frame
  write.csv(data3, file = paste0("cleaned_",new_name,".csv"), row.names = FALSE)
}
#sample files
fileNames_2 <- Sys.glob("cleaned_*")
samples_2 <- lapply(fileNames_2, function(x) {
  read.csv(x)
})

#Create Master_GP_list
GP_master_list <- lapply(fileNames_2, function(x) {
  read.csv(x)["GPs"]
})

# unlist to combine into a vector
GP_master_list <- data.frame(unlist(GP_master_list))
GP_master_list_cleaned <- unique(GP_master_list)
rownames(GP_master_list_cleaned) = NULL
colnames(GP_master_list_cleaned) <- "GPs"

# Create a new data frame to store the aligned data
aligned_df <- GP_master_list_cleaned

# Align data by GP master list, now renamed as aligned_df
for (df in samples_2) {
  aligned_df <- merge(aligned_df, df, all = TRUE)
}

# Filter data - the score of the GPs may be adjusted. Removing GPs with no scan
# time was used to reduce false discoveries.
aligned_df2 <- subset(aligned_df, aligned_df$Score > 300 & !is.na(aligned_df$Scan.Time))

duplicate_GPs <- subset(aligned_df2, duplicated(GPs))
unique_GPs <- subset(aligned_df2, !duplicated(GPs))

grouped <- aggregate(Scan.Time ~ GPs + Sample,data = aligned_df2, FUN = mean)

wider1 <- grouped %>%
pivot_wider(names_from = Sample, values_from = Scan.Time)

#GPs found in all samples
In_all <- na.omit(wider1)

# Export Results to Excel
write.xlsx(aligned_df2, "Byonic_align_results.xlsx", sheetName = "aligned_GPs")
write.xlsx(duplicate_GPs, "Byonic_align_results.xlsx", sheetName = "duplicates_or_isomers", append = TRUE)
write.xlsx(unique_GPs, "Byonic_align_results.xlsx", sheetName = "unique_IDs", append = TRUE)
write.xlsx(In_all, "Byonic_align_results.xlsx", sheetName = "found_in_all", append = TRUE)
