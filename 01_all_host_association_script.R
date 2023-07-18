# Creating a Merged Data Frame for Viral Contigs and MAGs

# Run this R script to create a csv file with all 4 sample's contigs and vMAGs combined into 1 dataset.

# Load packages
library(readr)
library(dplyr)

# Input data
jv119_contigs <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-119_874814_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv")
jv121_contigs <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-121_874818_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv")
jv132_contigs <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-132_874826_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv")
jv154_contigs <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-154_874822_Viral_Files/viral_hosts/viral_host_associations_filtered.tsv")
jv119_vMAG <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-119_874814_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv")
jv121_vMAG <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-121_874818_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv")
jv132_vMAG <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-132_874826_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv")
jv154_vMAG <- read_tsv("~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/jv-154_874822_Viral_Files/vMag_hosts/vMAG_host_associations_filtered.tsv")

# Create sample_name column for each data frame
jv119_contigs$sample_name <- "JV119"
jv119_vMAG$sample_name <- "JV119"
jv121_contigs$sample_name <- "JV121"
jv121_vMAG$sample_name <- "JV121"
jv132_contigs$sample_name <- "JV132"
jv132_vMAG$sample_name <- "JV132"
jv154_contigs$sample_name <- "JV154"
jv154_vMAG$sample_name <- "JV154"

# Create virus_type column for each data frame
jv119_contigs$virus_type <- "contig"
jv119_vMAG$element <- "vMAG"
jv121_contigs$element <- "contig"
jv121_vMAG$element <- "vMAG"
jv132_contigs$element <- "contig"
jv132_vMAG$element <- "vMAG"
jv154_contigs$element <- "contig"
jv154_vMAG$element <- "vMAG"

# Create sample_depth column for each data frame
jv119_contigs$sample_depth <- 400
jv119_vMAG$sample_depth <- 400
jv121_contigs$sample_depth <- 95
jv121_vMAG$sample_depth <- 95
jv132_contigs$sample_depth <- 80
jv132_vMAG$sample_depth <- 80
jv154_contigs$sample_depth <- 140
jv154_vMAG$sample_depth <- 140

# Create a list of all the data frames
df_list <- list(jv119_contigs, jv119_vMAG, jv121_contigs, jv121_vMAG, jv132_contigs, jv132_vMAG, jv154_contigs, jv154_vMAG)

# Create a list of the new columns so that all data frames have the same column names
new_names <- c("virus_name", "virus_length", "virus_read_count", "virus_read_depth", "virus_read_depth_in_host", "host_name","host_length", "host_read_count", "host_read_depth", "intra_read_count", "intra_linkage_density", "inter_read_count", "raw_inter_linkage_density", "raw_inter_vs_intra_ratio", "viral_copies_per_cell", "adjusted_inter_linkage_density", "adjusted_inter_vs_intra_ratio", "sample_name", "virus_type", "sample_depth")

# for loop to rename the columns in all of the data frames
for (i in seq_along(df_list)) {
  colnames(df_list[[i]]) <- new_names
}

# Use the dplyr package's bind_rows function to merge all data frames into one data frame
combo <- bind_rows(df_list[[1]], df_list[[2]], df_list[[3]], df_list[[4]], df_list[[5]], df_list[[6]], df_list[[7]], df_list[[8]]) # create a large data frame with all of the sample data frames combined

# Write the CSV file
#write.csv(combo,"~/Documents/Bigelow- starting Sept 2022/Virus Project/OMZ_MH_Analysis/Data/proximeta_viral_files/all_host_associations.csv" )
