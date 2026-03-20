##### Load package #####
library(RLDNe)
library(dplyr)

##### Give execution permissions to the Ne2-1L binary #####
ne_binary <- system.file("bin/linux/Ne2-1L", package = "RLDNe")
if (file.exists(ne_binary)) {
  Sys.chmod(ne_binary, mode = "0755")
}

##### Load arguments #####
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
apply_correction <- as.logical(args[12])
apply_harmo      <- as.logical(args[13])

next_args <- 14
if (apply_correction == TRUE) {
  nb_chrom <- as.numeric(args[next_args])
  next_args <- next_args + 1
}

if(apply_harmo == TRUE) {
  maf <- unlist(strsplit(args[next_args], "\\s+"))
  maf <- c(maf,"0+") #Systematic calculation for MAF 0+
}

##### Validate inputs #####
# Genepop files
if (length(list_gen) == 0) stop("gen_files is empty.", call. = FALSE)
if (length(list_gen) != length(list_names)) {
  stop("gen_files and gen_names must have the same number of entries ",
       "(got ", length(list_gen), " vs ", length(list_names), ").", call. = FALSE)
}

missing_gen <- list_gen[!file.exists(list_gen)]
if (length(missing_gen) > 0) {
  stop("Genepop input file(s) not found: ", paste(missing_gen, collapse = ", "), call. = FALSE)
}

#Indpop file
indpop <- tryCatch(
  read.table(indpop_path, header = FALSE, sep = "\t"),
  error = function(e) stop("Could not read indpop file: ", indpop_path,
                            "\n  ", conditionMessage(e), call. = FALSE)
)

if (ncol(indpop) < 2) {
  stop("indpop must have at least 2 tab-separated columns (Individual, Population). ",
       "Found ", ncol(indpop), " column(s).", call. = FALSE)
}

colnames(indpop) <- c("Ind", "Pop")


########################## RLDNe execution ##########################

#####################################################################
# Function : anonymise_genepop
# Description : Convert individuals name in number to avoid issue due to 
# the longer of individual's names.
######################################################################
anonymise_genepop <- function(gen_path) {
  ###### Change individuals' names to avoid reading problems for RLDNe #####
  lines <- readLines(gen_path)

  if (!file.exists(gen_path)) {
    warning(paste("File not found:", gen_path))
    next
  }

  # Find POP line
  pop_indices <- grep("^POP", lines, ignore.case = TRUE)
  if (length(pop_indices) != 1) {
    stop("Expected exactly one pop in genepop file")
  }

  # Copy lines (we will modify only individual names)
  new_lines <- lines

  # Counter for new individual IDs
  ind_counter <- 1

  # Store original individual names BEFORE renaming
  original_inds <- c()

  pop_index <- pop_indices[1]
  # Loop on individuals (lines after POP)
  for (k in (pop_index + 1):length(lines)) {

    current_line <- trimws(lines[k])

    # Skip empty lines if any
    if (nchar(current_line) == 0) next

    # Extract original individual name (before first comma)
    original_name <- trimws(sub(",.*", "", lines[k]))
    original_inds <- c(original_inds, original_name)

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
  list(path=n_gen_path, individuals=original_inds)
  #ajouter pour récupérer le nom de la pop --> stocker pour les sous-échantillons ??
}

######################################################################
# Function : extract_base_name
# Description : Extract the input real name. 
#####################################################################
extract_base_name <- function(input_name) {
  #Extract base name
  input_name <- basename(input_name)
  # Extract content from last parentheses if present
  if (grepl("\\([^)]+\\)\\s*$", input_name)) {
    input_name <- sub(".*\\(([^)]+)\\)\\s*$", "\\1", input_name)
  } else {
    input_name <- sub("\\.[^.]+$", "", input_name)
  }  
  return(input_name)
}

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
  #Remove extension
  base_params <- sub("\\.txt$", "", gen_name) 

  #Output file path
  results_file <- paste0("results_LDNe/", base_params, "_LDNe_results.txt")

  params_file_name <- paste0("LDNe_params_", base_params, ".txt")

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
results_individuals <- list()  # Store original individual names per result file
i_name <- 0
for (gen_path in list_gen) {
  #fonction 

  ##### Extract basename #####
  i_name <- i_name + 1
  gen_base <- extract_base_name(list_names[i_name])

  # Change the name of individuals to avoid issue due to the length of individual's name
  outputs <- anonymise_genepop(gen_path)
  n_gen_path <- outputs$path
  original_inds <- outputs$individuals

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
    results_individuals[[params_output$results_file]] <- original_inds

  ##### Run LDNe and return result file #####
  std_out <- RLDNe::run_LDNe(params_output$params_file)
}

############################# Extract RLDNe results ################################
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
# Population name is retrieved from indpop based on individuals
# present in the input genepop file (not from the filename)
###################################################################
extract_info_dataset <- function(result_path, indpop, individuals) {

  # Match individuals from the genepop to indpop to find the population
  matched_pops <- indpop$Pop[indpop$Ind %in% individuals]

  if (length(matched_pops) == 0) {
    stop(paste0(
      "No individuals from the genepop file were found in indpop for: ",
      result_path,
      "\nCheck that individual names match between the genepop and indpop files."
    ))
  }

  unique_pops <- unique(as.character(matched_pops))

  if (length(unique_pops) > 1) {
    warning(paste0(
      "Multiple populations found for file: ", result_path,
      " -> ", paste(unique_pops, collapse = ", "),
      ". Using the first one."
    ))
  }

  pop_in_filename <- unique_pops[1]

  #Extract basename (for dataset name and subsample)
  filename <- gsub("\\.txt", "",
                   basename(result_path), ignore.case = TRUE)

  #Extract name of the dataset
  dataset_name <- sub(paste0("_", pop_in_filename, ".*"), "", filename)
  # Fallback: if pop not in filename, use full basename
  if (dataset_name == filename) {
    dataset_name <- filename
  }

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
  lines <- readLines(text_file)

  # Check result from RLDNe before extract data
  if (any(grepl("Fatal error", lines, ignore.case = TRUE))) {
    message("File ignored due to Fatal error: ", text_file)
    next
  }

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
  info <- extract_info_dataset(text_file, indpop, results_individuals[[text_file]])
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

######################## Pseudoreplication correction ######################
if (apply_correction == TRUE) {
  ldne_results <- ldne_results %>%
    mutate(
      NeLD = as.numeric(NeLD),
      NeLD = ifelse(is.na(NeLD), 999999, NeLD),
      NeLD_corrected = ifelse(
        NeLD == 999999,
        999999,
        NeLD / (0.098 + 0.219 * log(nb_chrom))
      )
    )
}

write.table(ldne_results,
            file = "summary_file/LDNe_results.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

######################## Harmonic mean ###########################
if (apply_harmo == TRUE) {
  harmonic_mean <- function(x) {
    n <- length(x)
    x[is.na(x)] <- 999999 #NA means Infinite
    return(n / sum(1 / x))
  }

  #Detect which colums are present in the input file
  available_cols <- colnames(ldne_results)

  numeric_cols <- c("NeLD", "JK_CI_down",
                  "JK_CI_up",
                  "Overall_LD_r2",
                  "Expected_LD_r2",
                  "NeLD_corrected")
  cols_to_use <- numeric_cols[numeric_cols %in% available_cols]

  #Final table
  ne_estim_all <- data.frame()

  for (i in maf){
    ldne_harm <- ldne_results
    for (col in cols_to_use) {
      ldne_harm[[col]] <- as.numeric(ldne_harm[[col]])
    }

    # Count the number of subsets for each Pop and Marker_type combination
    subset_counts <- ldne_harm %>%
      filter(MAF == i) %>%
      group_by(Dataset, Pop, Marker_type) %>%
      summarize(n_subsets = n(), .groups = "drop")
  
    #Build the summarize file
    ne_estim <- ldne_harm %>%
      filter(MAF == i) %>%
      group_by(Dataset, Pop, Marker_type) %>%
      summarize(across(all_of(cols_to_use), harmonic_mean))
  
    ne_estim$MAF <- i
  
    ne_estim <- ne_estim %>%
      left_join(subset_counts, by = c("Dataset", "Pop", "Marker_type")) %>%
      mutate(Subset = paste0("H_mean_btw_", n_subsets, "sub")) %>%
      select(-n_subsets)  # Remove the temporary count column
      relocate(Subset, .after = "Marker_type")

      ne_estim_all <- bind_rows(ne_estim_all, ne_estim)
  }

  output_file <- paste0("summary_file/LDNe_harmonic_mean.txt")

  # save file
  write.table(ne_estim_all, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
