#Load packages
library(dplyr)

#Load arguments
args <- commandArgs(trailingOnly = TRUE)

table_results_path <- args[1]
results <- read.table(table_results_path, header = TRUE, sep = "\t")
maf <- unlist(strsplit(args[2], "\\s+"))
maf <- c(maf,"0+") #Systematic calculation for MAF 0+

#Create output directory
dir.create("summary_file", recursive = TRUE)

####################################################################
# Function : harmonic_mean
# Description Function to calculate harmonic mean
####################################################################
harmonic_mean <- function(x) {
  n <- length(x)
  x[is.na(x)] <- 999999 #NA means Infinite
  return(n / sum(1 / x))
}

#############################
# Calculate the harmonic mean
#############################

#Detect which colums are present in the input file
available_cols <- colnames(results)

numeric_cols <- c("NeLD", "JK_CI_down",
                  "JK_CI_up",
                  "Overall_LD_r2",
                  "Expected_LD_r2",
                  "NeLD_corrected")
cols_to_use <- numeric_cols[numeric_cols %in% available_cols]

#Final table
ne_estim_all <- data.frame()

for (i in maf){
  results2 <- results
  for (col in cols_to_use) {
    results2[[col]] <- as.numeric(results2[[col]])
  }

  # Count the number of subsets for each Pop and Marker_type combination
  subset_counts <- results2 %>%
    filter(MAF == i) %>%
    group_by(Dataset, Pop, Marker_type) %>%
    summarize(n_subsets = n(), .groups = "drop")
  
  #Build the summarize file
  ne_estim <- results2 %>%
    filter(MAF == i) %>%
    group_by(Dataset, Pop, Marker_type) %>%
    summarize(across(all_of(cols_to_use), harmonic_mean))
  
  ne_estim$MAF <- i
  
  ne_estim <- ne_estim %>%
    left_join(subset_counts, by = c("Dataset", "Pop", "Marker_type")) %>%
    mutate(Subset = paste0("H_mean_btw_", n_subsets, "sub")) %>%
    select(-n_subsets)  # Remove the temporary count column

  ne_estim_all <- bind_rows(ne_estim_all, ne_estim)
}

output_file <- paste0("summary_file/LDNe_harmonic_mean.txt")

# save file
write.table(ne_estim_all, output_file, sep = "\t", row.names = FALSE, quote = FALSE)