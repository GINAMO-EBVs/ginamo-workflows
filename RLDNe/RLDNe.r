#Load package
library(RLDNe)
library(dplyr)

# Give execution permissions to the Ne2-1L binary
ne_binary <- system.file("bin/linux/Ne2-1L", package = "RLDNe")
if (file.exists(ne_binary)) {
  Sys.chmod(ne_binary, mode = "0755")
}

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
list_gen <- unlist(strsplit(args[1], ",\\s*"))
list_names <- unlist(strsplit(args[2], ",\\s*"))
stopifnot(length(list_gen) == length(list_names))

indpop_path <- args[3]
indpop <- read.table(indpop_path, header = FALSE, sep = "\t")
colnames(indpop) <- c("Ind", "Pop")
marker_type <- args[4]
n_critical_values <- args[5]
critical_freqs <- args[6]
tabular_output <- 0
confidence_intervals <- args[7]
mating_system <- args[8]
max_individuals <- args[9]
pop_range <- args[10]
loc_range <- args[11]
ld_method <- 1

#Transform list_gen into a vector to implement across all files
list_gen <- unlist(strsplit(args[1], ",\\s*"))

# Create output directory
dir.create("results_LDNe", recursive = TRUE)
dir.create("genfiles", recursive = TRUE)

########################## RLDNe execution ##########################

#######################################################################
# Function : params_file
# Description : Creation of the parameter file required for NeEstimator
#######################################################################
params_file <- function(gen_path,
                        gen_name,
                        ld_method,
                        n_critical_values,
                        critical_freqs,
                        tabular_output,
                        confidence_intervals,
                        mating_system,
                        max_individuals,
                        pop_range,
                        loc_range) {
  #Extract basename
  base_name <- basename(gen_name)
  base_name <- sub("\\.txt$", "", base_name) #remove extension

  #Output file path
  results_file <- paste0("results_LDNe/", base_name, "_LDNe_results.txt")

  params_file_name <- paste0("LDNe_params_", base_name, ".txt")

  lines <- c(
    paste0(ld_method, "\t* LD Method"),
    paste0(n_critical_values, "\t* number of critical values"),
    paste0(critical_freqs, "\t* critical allele frequency values"),
    paste0(tabular_output, "\t* tabular output"),
    paste0(confidence_intervals, "\t* confidence intervals"),
    paste0(mating_system, "\t* 0: Random mating, 1: Monogamy (LD method)"),
    paste0(max_individuals,
           "\t* max individual to be processed per pop, 0 for no limit"),
    paste0(pop_range,
           "\t* Pop. range to run, given in pair. No limit if the first = 0"),
    paste0(loc_range,
           "\t* Loc. ranges to run, given in pairs. No limit if the 1st = 0"),
    paste0(results_file, "\t* output file name"),
    paste0(gen_path, "\t* input file")
  )

  writeLines(lines, params_file_name)
  return(list(params_file = params_file_name, results_file = results_file))
}

##########################
# LDNe : Main execution
##########################
results_path <- c()
i_name <- 0
for (gen_path in list_gen) {
  if (!file.exists(gen_path)) {
    warning(paste("File not found:", gen_path))
    next
  }

  ###### Change individuals' names to avoid reading problems for RLDNe #####
  lines <- readLines(gen_path)

  # Find POP line
  pop_indices <- grep("^POP", lines, ignore.case = TRUE)
  if (length(pop_indices) != 1) {
    stop("Not a single pop in this dataset, check please")
  }

  pop_index <- pop_indices[1]

  # Copy lines (we will modify only individual names)
  new_lines <- lines

  # Counter for new individual IDs
  ind_counter <- 1

  # Loop on individuals (lines after POP)
  for (k in (pop_index + 1):length(lines)) {

    current_line <- trimws(lines[k])

    # Skip empty lines if any
    if (nchar(current_line) == 0) next

    # Replace only what is before the first comma
    new_lines[k] <- sub(
      "^[^,]+",
      ind_counter,
      lines[k]
    )

    ind_counter <- ind_counter + 1
  }

  # Build output file name
  name_file <- sub("^.*/(.*)\\.[^.]+$", "\\1", gen_path)
  n_gen_path <- paste0("genfiles/", name_file, ".txt")

  # Save file
  writeLines(new_lines, n_gen_path)

  if (!file.exists(n_gen_path)) {
    warning(paste("File not found:", n_gen_path))
    next
  }

  ##### Extract basename #####
  i_name <- i_name + 1

  #Extract base name
  # Extract content from last parentheses if present
  if (grepl("\\([^)]+\\)\\s*$", list_names[i_name])) {
    # Extract content between last parentheses
    gen_base <- sub(".*\\(([^)]+)\\)\\s*$", "\\1", list_names[i_name])
  } else {
    # No parentheses, use the name without extension
    gen_base <- list_names[i_name]
  }

  ###### Create params file #####
  params_output <- params_file(n_gen_path,
                      gen_base,
                        ld_method,
                        n_critical_values,
                        critical_freqs,
                        tabular_output,
                        confidence_intervals,
                        mating_system,
                        max_individuals,
                        pop_range,
                        loc_range)

    ##### Extract params outputs #####
    results_path <- c(results_path, params_output$results_file) 
    
  ##### Run LDNe and return result file #####
  std_out <- RLDNe::run_LDNe(params_output$params_file)
}

############################# Extract RLDNe results ################################

#Create output directory
dir.create("summary_file", recursive = TRUE)

#Create output file
ldne_results <- as.data.frame(matrix(ncol = 11, nrow = 0))
colnames(ldne_results) <- c("Dataset",
                            "Marker_type",
                            "Pop",
                            "Subset",
                            "N_loci",
                            "MAF",
                            "NeLD",
                            "JK_CI_down",
                            "JK_CI_up",
                            "Overall_LD_r2",
                            "Expected_LD_r2")

###################################################################
# Function : extract_info_dataset
# Description : Function to extract information about the dataset
# (pop, dataset name and subsample number)
###################################################################
#Important : population name need to be on file name
extract_info_dataset <- function(result_path, indpop) {
  #Extract populations name
  pops <- as.character(unique(indpop[,2]))

  #Extract basename
  filename <- gsub("\\.txt", "",
                   basename(result_path), ignore.case = TRUE)

  #Find which pop is present on dataset name
  pop_in_filename <- pops[sapply(pops, function(p) grepl(p, filename))]

  # If no pop is found → informative error
  if (length(pop_in_filename) == 0) {
    stop("No pop found in the file name. 
    Check the indpop file or the dataset name.")
  }

  #Extract name of the dataset
  dataset_name <- sub(paste0("_", pop_in_filename, ".*"), "", filename)

  #Extract subsample number
  subsample_match <- regmatches(filename, regexpr("subsample_[0-9]+", filename))

  if (length(subsample_match) == 1) {
    subsample_number <- sub("subsample_", "", subsample_match)
  } else {
    subsample_number <- NA
  }

  return(list(
    dataset_name = dataset_name,
    pop_in_filename  = pop_in_filename,
    subsample = subsample_number
  ))
}

###################################################################
# Function : extract_values
# Description : allows you to extract values from RLDNe result file
###################################################################

extract_values <- function(line) {
  values <- unlist(strsplit(line, "\\s+"))
  values <- values[values != ""]
}

###################################################################
# Extract value execution
##################################################################
n_critical_values=as.numeric(n_critical_values)
n_critical_values=n_critical_values+1
for (text_file in results_path){
  # Check result from RLDNe before extract data
  if (any(grepl("Fatal error", lines, ignore.case = TRUE))) {
    message("File ignored due to Fatal error: ", text_file)
    next
  }

  lines <- readLines(text_file)

  #Lines with interesting data
  start_index <- grep("Lowest Allele Frequency Used", lines)
  end_index <- grep("Ending time:", lines)

  #Get interesting lines
  table_lines <- lines[(start_index + 4):(end_index - 4)]

  #Remove lines not interesting
  table_lines <- table_lines[table_lines != ""]

  #Extract values
  estimated_ne <- extract_values(table_lines[5])
  jk_ci_down <- extract_values(table_lines[9])
  jk_ci_up <- extract_values(table_lines[10])
  overall_ld_r2 <- extract_values(table_lines[3])
  expected_ld_r2 <- extract_values(table_lines[4])
  crit_freqs <- unlist(strsplit(args[6], "\\s+"))
  crit_freqs <- c(crit_freqs, "0+")

  #Get information about loci used to estimate LDNe
  snp_line_index <- grep("Number of Loci", lines)
  snp_line <- lines[snp_line_index]
  snp_used <- extract_values(snp_line[1])[5]

  #Get dataset and population name
  info <- extract_info_dataset(text_file, indpop)
  dataset <- info$dataset_name
  population <- info$pop_in_filename
  subsample <- info$subsample

  #Store information in a tab
  ldne_results <- rbind(ldne_results,
                        data.frame(Dataset = rep(dataset, n_critical_values),
                                   Pop = rep(population, n_critical_values),
                                   Marker_type = rep(marker_type, n_critical_values),
                                   Subset = rep(subsample, n_critical_values),
                                   loci = rep(snp_used, n_critical_values),
                                   MAF = crit_freqs,
                                   NeLD = tail(estimated_ne, n_critical_values),
                                   JK_CI_down = tail(jk_ci_down,
                                                     n_critical_values),
                                   JK_CI_up = tail(jk_ci_up, n_critical_values),
                                   Overall_LD_r2 = tail(overall_ld_r2,
                                                        n_critical_values),
                                   Expected_LD_r2 = tail(expected_ld_r2,
                                                         n_critical_values)))

}

ldne_results <- ldne_results %>%
  mutate(across(c("NeLD",
                  "JK_CI_down",
                  "JK_CI_up",
                  "Overall_LD_r2",
                  "Expected_LD_r2"), as.numeric),
    JK_CI_up = ifelse(is.na(JK_CI_up), 999999, JK_CI_up),
    NeLD     = ifelse(is.na(NeLD), 999999, NeLD)
  )


write.table(ldne_results,
            file = "summary_file/LDNe_results.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
