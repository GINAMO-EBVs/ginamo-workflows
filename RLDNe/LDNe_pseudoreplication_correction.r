#Load packages
library(dplyr)

#Load arguments
args <- commandArgs(trailingOnly = TRUE)

table_results_path <- args[1]
results <- read.table(table_results_path, header = TRUE, sep = "\t")
nb_chrom <- as.numeric(args[2])

#Create output directory
dir.create("summary_file", recursive = TRUE)

# Correction
results <- results %>%
  mutate(
    NeLD = as.numeric(NeLD),
    NeLD = ifelse(is.na(NeLD), 999999, NeLD),
    NeLD_corrected = ifelse(
      NeLD == 999999,
      999999,
      NeLD / (0.098 + 0.219 * log(nb_chrom))
    )
  )

output_file <- paste0("summary_file/LDNe_pseudoreplication_corrected.txt")

# save file
write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)