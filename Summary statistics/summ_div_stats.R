#Load packages
library(dplyr)
library(ade4)
library(adegenet)
library(vcfR)
library(hierfstat)
library(mmod)

#Load arguments
args <- commandArgs(trailingOnly = TRUE)

marker_type <- args[1]
input_path  <- args[2]
input_name  <- args[3]

# Load indpop only for SNP
if (marker_type == "snp") {
  indpop_path <- args[4]
  indpop <- read.table(indpop_path, header = FALSE, sep = "\t")
  colnames(indpop) <- c("Ind", "Pop")
}

#Create output file
summary_div_stats <- as.data.frame(matrix(ncol = 11, nrow = 0))
colnames(summary_div_stats) <- c("Dataset",
                                 "Pop",
                                 "Marker_type",
                                 "n_ind",
                                 "Hobs",
                                 "Hexp",
                                 "Ar",
                                 "Fis",
                                 "average_pairwise_Fst",
                                 "average_Gst_Nei",
                                 "average_Jost_D")

#############################################################################
# Function : load_snp_data
# Description : Loads SNP data from VCF file and converts to genind object
# with population assignments from indpop file.
#############################################################################

load_snp_data <- function(vcf_path, indpop) {
  #Load VCF file

  vcf_file <- vcfR::read.vcfR(vcf_path)

  #Conversion to genind object
  gen_file <- vcfR::vcfR2genind(vcf_file)

  # Add individuals and their population
  ind_list <- data.frame(Ind = adegenet::indNames(gen_file))
  colnames(ind_list) <- "Ind"

  # Check for individuals in VCF but not in indpop
  missing_ind <- setdiff(ind_list$Ind, indpop$Ind)
  if (length(missing_ind) > 0) {
    stop(sprintf("%d individu(s) du VCF absent(s) d'indpop : %s", 
                 length(missing_ind), 
                 paste(head(missing_ind, 5), collapse = ", ")))
  }

  ind_df <- ind_list %>%
    mutate(order = row_number()) %>%
    left_join(indpop, by = "Ind") %>%
    arrange(order) %>%
    select(-order) #make sure same order is kept

    adegenet::pop(gen_file) <- as.factor(ind_df$Pop) #add pop on genind object

    return(gen_file)

}

#############################################################################
# Function : load_ssr_data
# Description : Loads SSR/microsatellite data from tabular file and
# converts to genind object
#############################################################################

load_ssr_data <- function(input_path) {
  ssr_raw <- read.csv(input_path, sep="\t")

  geno_cols <- setdiff(colnames(ssr_raw), c("Ind", "Pop"))

  gen_file <- df2genind(ssr_raw[geno_cols],
                      pop = ssr_raw$Pop,
                      ind.names = ssr_raw$Ind,
                      ploidy = 2,
                      NA.char = "0",
                      sep = "/")

  
  return(gen_file)
}

##############################################################################
# Function : compute_div_stats
# Description : Core function that computes diversity statistics from a
# genind object. Works for both SNP and SSR data.
# Calculates: Hobs, Hexp, Fis and allelic richness per population.
##############################################################################

compute_div_stats <- function(gen_file, dataset_name, marker_type) {

  pops <- seppop(gen_file, drop = TRUE)

  # basic.stats calculation by pop
  bs_list <- lapply(pops, function(ls) {
    hf <- hierfstat::genind2hierfstat(ls)
    hierfstat::basic.stats(hf)
  })

  # Hobs, Hexp
  Hobs <- sapply(bs_list, function(bs) bs$overall["Ho"])
  Hexp <- sapply(bs_list, function(bs) bs$overall["Ht"])

  #Fis
  Fis <- t(sapply(seppop(gen_file), function(ls) basic.stats(ls)$perloc$Fis))
  Fis_pop <- rowMeans(as.matrix(Fis),na.rm=TRUE)

  # Allelic richness
  Richness <- hierfstat::allelic.richness(gen_file, diploid = TRUE)
  Richness_mean <- colMeans(as.matrix(Richness$Ar), na.rm = TRUE)

  Richness_mean <- Richness_mean[names(Richness_mean) != "dumpop"]
  pop_names <- levels(pop(gen_file))
  names(Richness_mean) <- pop_names


  # Nombre d'individus
  n_ind_pop <- table(pop(gen_file))

  div_stats_dataset <- data.frame(
    Dataset = dataset_name,
    Pop = names(n_ind_pop),
    Marker_type = marker_type,
    n_ind = as.vector(n_ind_pop),
    Hobs = as.vector(Hobs),
    Hexp = as.vector(Hexp),
    Ar   = Richness_mean,
    Fis  = as.vector(Fis_pop),
    average_pairwise_Fst = NA,
    average_Gst_Nei = NA,
    average_Jost_D = NA
  )

  return(div_stats_dataset)
}

########################################################################
# Function : pairwise_values_Fst_DJost
# Description : 
########################################################################
pairwise_values_Fst_DJost <- function(gen_path) {
  #Pairwise values
  fst <- genet.dist(gen_path, diploid = TRUE, method = "WC84")

  #Pairwise value from mmod package
  gst_pr_nei <- pairwise_Gst_Nei(gen_path)
  jost_D <- pairwise_D(gen_path)

  matrices_list <- list(
    fst = fst,
    gst_pr_nei = gst_pr_nei,
    jost_D = jost_D
  )
  
  # Remove temporary matrices
  rm(fst, gst_pr_nei, jost_D)
  
  return(matrices_list)
}


########################################################################
# Function : average_pairwise_by_pop
# Description : 
########################################################################

average_pairwise_by_pop <- function(matrices_list, dataset_name) {
  stats_diff <- c("fst", "gst_pr_nei", "jost_D")

  # Get population labels from first matrix
  labels <- attr(matrices_list[["fst"]], "Labels")
  
  # Initialize lists to store results
  pops_list <- c()
  fst_list <- c()
  gst_nei_list <- c()
  jost_d_list <- c()
  
  for (pop in labels) {
    mean_values <- c()
    
    for (measure in stats_diff) {
      values <- as.matrix(matrices_list[[measure]])
      all_pairwise_1pop <- values[pop, ]
      all_pairwise_1pop <- all_pairwise_1pop[names(all_pairwise_1pop) != pop] #to remove pop itself
      mean_values <- c(mean_values, mean(all_pairwise_1pop, na.rm = TRUE))
    }
    
    pops_list <- c(pops_list, pop)
    fst_list <- c(fst_list, mean_values[1])
    gst_nei_list <- c(gst_nei_list, mean_values[2])
    jost_d_list <- c(jost_d_list, mean_values[3])
  }
  
  result <- data.frame(
    Dataset = dataset_name,
    Pop = pops_list,
    average_pairwise_Fst = fst_list,
    average_Gst_Nei = gst_nei_list,
    average_Jost_D = jost_d_list,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

##############################################################
# Function : save_matrices
# Description : Save pairwise matrices to files
###############################################################

save_matrices <- function(matrices_list, dataset_name) {
  stats_diff <- c("fst", "gst_pr_nei", "jost_D")
  
  for (measure in stats_diff) {
    matrix_data <- as.matrix(matrices_list[[measure]])
    labels <- attr(matrices_list[[measure]], "Labels")
    
    # Add row and column names
    rownames(matrix_data) <- labels
    colnames(matrix_data) <- labels
    
    # Save to file
    filename <- paste0("matrices_output/", dataset_name, "_", measure, "_matrix.txt")
    write.table(matrix_data, 
                file = filename, 
                quote = FALSE, 
                sep = "\t",
                row.names = TRUE,
                col.names = NA)  # col.names = NA to align header with data
  }
}

#######################################
# Main execution 
#######################################

# Extract base name (remove extension, handle parentheses from Galaxy labels)
if (grepl("\\([^)]+\\)\\s*$", input_name)) {
  input_name <- sub(".*\\(([^)]+)\\)\\s*$", "\\1", input_name)
} else {
  input_name <- sub("\\.[^.]+$", "", input_name)
}  

# Load data according to marker type
if (marker_type == "snp") {
  gen_file <- load_snp_data(input_path, indpop)
} else if (marker_type == "ssr") {
  gen_file <- load_ssr_data(input_path)
}
  
######### By pop ############
# Compute diversity statistics
result <- compute_div_stats(gen_file, input_name, marker_type)

#Pairwise value (Fst, DJost, Gst_Nei)
pairwise_matrices <- pairwise_values_Fst_DJost(gen_file)

#Average pairwise for a population
mean_pairwise <- average_pairwise_by_pop(pairwise_matrices, dataset_name= input_name)

# Merge mean pairwise values with diversity stats
result <- result %>%
  left_join(mean_pairwise, by = "Pop", suffix = c("", "_new"))
  
# Update the mean columns with the new values
result$average_pairwise_Fst <- result$average_pairwise_Fst_new
result$average_Gst_Nei <- result$average_Gst_Nei_new
result$average_Jost_D <- result$average_Jost_D_new
  
# Remove duplicate columns
result <- result %>%
  select(Dataset, Pop, Marker_type, n_ind, Hobs, Hexp, Ar, Fis, average_pairwise_Fst, average_Gst_Nei, average_Jost_D)
  
# Add results to the global dataframe
if (!is.null(result) && nrow(result) > 0) {
  summary_div_stats <- rbind(summary_div_stats, result)
}

# Save matrices to files
save_matrices(pairwise_matrices, input_name)

############ Total (all populations combined) #################
# Reload data with "tot" population
if (marker_type == "snp") {
  #Create a modified indpop with "tot" for all individuals
  indpop_tot <- indpop
  indpop_tot$Pop <- "tot"
  gen_file_tot <- load_snp_data(input_path, indpop_tot)
} else if (marker_type == "ssr") {
  gen_file_tot <- load_ssr_data(input_path)
  # For SSR, manually set all populations to "tot"
  adegenet::pop(gen_file_tot) <- as.factor(rep("tot", nInd(gen_file_tot)))
}

# Compute diversity statistics for "tot"
result_tot <- compute_div_stats(gen_file_tot, input_name, marker_type)

#Mean pairwise values can't be estimate because only 1 population
result_tot$average_pairwise_Fst <- NA
result_tot$average_Gst_Nei <- NA
result_tot$average_Jost_D <- NA

# Add "tot results" to the global dataframe
if (!is.null(result_tot) && nrow(result_tot) > 0) {
  summary_div_stats <- rbind(summary_div_stats, result_tot)
}

########### MEAN (average across populations, excluding "tot" rows) #########
summary_means <- summary_div_stats %>%
  filter(Pop != "tot") %>%  # Exclude lines "tot"
  group_by(Dataset, Marker_type) %>%
  summarise(
    Pop = "mean",
    n_ind = sum(n_ind, na.rm = TRUE),
    Hobs = mean(Hobs, na.rm = TRUE),
    Hexp = mean(Hexp, na.rm = TRUE),
    Ar = mean(Ar, na.rm = TRUE),
    Fis = mean(Fis, na.rm = TRUE),
    average_pairwise_Fst = mean(average_pairwise_Fst, na.rm = TRUE),
    average_Gst_Nei = mean(average_Gst_Nei, na.rm = TRUE),
    average_Jost_D = mean(average_Jost_D, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(Dataset, Pop, Marker_type, n_ind, Hobs, Hexp, Ar, Fis, average_pairwise_Fst, average_Gst_Nei, average_Jost_D)

# Add to the final table
summary_div_stats <- rbind(summary_div_stats, summary_means)

#Round to 5 decimal
cols_to_round <- c("Hobs", "Hexp", "Ar", "Fis", "average_pairwise_Fst", "average_Gst_Nei", "average_Jost_D")
summary_div_stats[cols_to_round] <- lapply(summary_div_stats[cols_to_round], round, 5)


write.table(summary_div_stats,
              file = "summary_file/summary_div_stats.txt",
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")

