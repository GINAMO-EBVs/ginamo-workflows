# Load packages
library(adegenet)
library(poppr)
library(hierfstat)
 
# Parse named arguments  --flag value
args <- commandArgs(trailingOnly = TRUE)
 
get_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) stop("Missing argument: ", flag)
  return(args[idx + 1])
}
 
input_path     <- get_arg(args, "--input")
input_name     <- get_arg(args, "--name")
filter_order   <- get_arg(args, "--order")
ind_md         <- as.numeric(get_arg(args, "--ind-md"))
ssr_md         <- as.numeric(get_arg(args, "--ssr-md"))
low_quar       <- as.numeric(get_arg(args, "--low-quar"))
high_quar      <- as.numeric(get_arg(args, "--high-quar"))
multi          <- as.numeric(get_arg(args, "--multi"))
threshold_pval <- as.numeric(get_arg(args, "--threshold-pval"))
 
# Clean input_name (handle Galaxy element_identifier with parentheses)
if (grepl("\\([^)]+\\)\\s*$", input_name)) {
  input_name <- sub(".*\\(([^)]+)\\)\\s*$", "\\1", input_name)
} else {
  input_name <- sub("\\.[^.]+$", "", input_name)
}
 
# Parse filter order from comma-separated string (e.g. "A,B,C,D")
# Map letters to internal filter names, preserving user-defined order
filter_map <- list(
  A = "ind_md_ssr",
  B = "ssr_md_filtering",
  C = "null_allele",
  D = "LD_filtering"
)
 
if (filter_order == "") {
  stop("No filter selected. Please select at least one filter to apply.")
}

selected_letters <- trimws(strsplit(filter_order, ",")[[1]])
filters <- unlist(filter_map[selected_letters], use.names = FALSE)
filters <- filters[!is.na(filters)]
 
message("Filters to apply (in order): ", paste(filters, collapse = " -> "))
 

# Load data
data <- read.csv(input_path, head = TRUE, sep = "\t")

# Genind conversion
data_genind <- data[, !colnames(data) %in% c("Ind", "Pop")]
genind_file <- df2genind(
  data_genind,
  sep = "/",
  NA.char = "0/0",
  ploidy = 2,
  pop = data$Pop,
  ind.names = data$Ind
)


#####################################################################################
# Function : ind_md_filtering
# Description : Remove individuals with more missing data than the defined threshold
######################################################################################

ind_md_filtering <- function(genind_file, ind_md){
  genind_ind_md <- missingno(genind_file, type = "loci", cutoff = ind_md)

  return(genind_ind_md)
}

####################################################################################
# Function : ssr_md_filtering
# Description : Remove SSR with more missing data than the defined threshold
####################################################################################

ssr_md_filtering <- function(genind_file, ssr_md){
  genind_ssr_md <- missingno(genind_file, type = "geno", cutoff = ssr_md)

  return(genind_ssr_md)
}

######################################################################################################################
# Function : null_allele
# Description : Detection and removal of null or paralogue alleles using the distribution of Fis. 
#               Outliers are considered to be probable null alleles and are removed. 
#               For the detection of outliers, the threshold is estimate using the IQR (interquartile range) method.
#####################################################################################################################

null_allele <- function(genind_file, low_quar, high_quar, multi) {
  # Fis calculation
  basic_stats <- basic.stats(genind_file)

  # Extraction Fis per loci
  fis_per_locus <- basic_stats$perloc$Fis

  # Boxplot visualisation 
  png_name=paste0("figure_output/", input_name, "_fis_null_allele.png")
  png(png_name, width = 800, height = 600)
  boxplot(fis_per_locus,
          main = paste0("Fis distribution per locus: ", input_name),
          ylab = "Fis",
          col  = "steelblue")
  dev.off()

  # Statistical identification of outliers (IQR method)
  quartiles <- quantile(fis_per_locus, probs=c(low_quar, high_quar), na.rm = TRUE)
  IQR <- diff(quartiles)

  low_threshold <- quartiles[1] - multi * IQR
  high_threshold <- quartiles[2] + multi * IQR

  # Identify the names of outlier loci
  outliers_loci <- names(fis_per_locus[fis_per_locus < low_threshold | fis_per_locus > high_threshold])
  print(paste("Locus outliers identified :", paste(outliers_loci, collapse=", ")))

  # Only keep the loci that are NOT in the list of outliers.
  locus_to_keep <- setdiff(locNames(genind_file), outliers_loci)
  genind_without_outliers <- genind_file[loc = locus_to_keep]

  return(genind_without_outliers)
}

####################################################################
# Function : genind_to_genepop
# Description : Convert genind into a genepop file, compatible with 
#               Genepop package.
#####################################################################

genind_to_genepop <- function(genind_file, file) {
  
  tab <- genind_file@tab
  loci <- genind_file@loc.fac
  loc_names <- locNames(genind_file)
  pops <- pop(genind_file)
  
  # Automatic detection of the number of digits describing each allle
  all_alleles <- unlist(alleles(genind_file))
  digit_num <- max(nchar(all_alleles), na.rm = TRUE)
  message("Number of digits describing each allele : ", digit_num)
  
  lines_out <- c("SSR data converted from genind")
  lines_out <- c(lines_out, loc_names)
  
  for (p in levels(pops)) {
    lines_out <- c(lines_out, "POP")
    inds <- which(pops == p)
    
    for (i in inds) {
      geno_parts <- c()
      
      for (loc in loc_names) {
        alleles <- tab[i, loci == loc]
        allele_names <- gsub(".*\\.", "", names(alleles))
        called <- allele_names[!is.na(alleles) & alleles > 0]
        
        if (length(called) == 0) {
          geno_parts <- c(geno_parts, strrep("0", digit_num * 2))
        } else if (length(called) == 1) {
          a <- sprintf(paste0("%0", digit_num, "d"), as.integer(called[1]))
          geno_parts <- c(geno_parts, paste0(a, a))
        } else {
          a1 <- sprintf(paste0("%0", digit_num, "d"), as.integer(called[1]))
          a2 <- sprintf(paste0("%0", digit_num, "d"), as.integer(called[2]))
          geno_parts <- c(geno_parts, paste0(a1, a2))
        }
      }
      
      ind_name <- rownames(tab)[i]
      lines_out <- c(lines_out, paste0(ind_name, " ,  ", paste(geno_parts, collapse = " ")))
    }
  }
  
  writeLines(lines_out, file)
}


#################################################################################################
# Function : LD_filtering
# Description : Performs the exact test for each pair of loci for each population using Genepop 
#               Deletes the locus involved  in the most pairs with significant p-values.
#################################################################################################

LD_filtering <- function(genepop_file_path, threshold_pval, genind_file) {
 
  output_ld <- "LD_test.txt"
 
  # Write settings file for the Genepop binary
  # Option 2 = Linkage disequilibrium exact test
  settings_file <- "genepop_LD_settings.txt"
  writeLines(c(
    paste0("InputFile=", genepop_file_path),
    paste0("OutputFile=", output_ld),
    "MenuOptions=2.1",
    "Dememorization=10000",
    "BatchNumber=100",
    "BatchSize=5000"
  ), settings_file)
 
  # Call the Genepop binary
  ret <- system2("Genepop", args = settings_file, stdout = TRUE, stderr = TRUE)
  message(paste(ret, collapse = "\n"))
 
  if (!file.exists(output_ld)) {
    stop("Genepop binary did not produce the output file. Check the log above.")
  }
 
  # Read output file from Genepop binary
  # Output format is identical to what test_LD() from the R package produces
  lines <- readLines(output_ld)
  lines <- trimws(lines, which = "right")
 
  # Extraction of the pairs with significant p-values
  # Extract only the lines containing the results for locus pairs
  global_lines <- lines[grepl("&", lines)]
 
  global_table <- do.call(rbind, lapply(global_lines, function(l) {
    parts <- strsplit(trimws(l), "\\s+")[[1]]
    # Extract p-values from the 6th position and remove special characters
    pval_raw <- gsub("[><]", "", parts[6])
    # Safety check: if the 6th part is empty, try to get the P-value from the 5th part
    if (is.na(pval_raw)) pval_raw <- gsub("[><]", "", parts[5])
    data.frame(
      Locus1  = parts[1],
      Locus2  = parts[3],
      P_value = as.numeric(pval_raw),
      stringsAsFactors = FALSE
    )
  }))
 
  # Filter the table to keep only pairs with significant P-values
  paires_sig <- global_table[global_table$P_value < threshold_pval, ]
 
  # Greedy strategy : Deletes the locus involved in the most pairs with significant p-values.
  loci_to_remove <- c()
  paires_remaining <- paires_sig
 
  while (nrow(paires_remaining) > 0) {
  
    # Count how many times each locus appears in the remaining pairs
    tous_loci <- c(paires_remaining$Locus1, paires_remaining$Locus2)
    compte <- sort(table(tous_loci), decreasing = TRUE)
  
    # Remove the most frequent locus
    locus_to_remove <- names(compte)[1]
    loci_to_remove <- c(loci_to_remove, locus_to_remove)
  
    # Remove all pairs involving this locus
    paires_remaining <- paires_remaining[
      paires_remaining$Locus1 != locus_to_remove & 
        paires_remaining$Locus2 != locus_to_remove, 
    ]
  }
 
  message("Deleted loci : ", paste(loci_to_remove, collapse = ", "))
 
  # Delete loci in the genind
  loci_to_keep <- setdiff(locNames(genind_file), loci_to_remove)
  genind_nLD <- genind_file[loc = loci_to_keep]
 
  return(genind_nLD)
 
}

########################
# Main execution
# Conversion et test LD sur le genind final (après filtrage missing data et allèles nuls)
##########################

# Summary table
summary_rows <- list()
 
# Filter labels and parameter strings
filter_labels <- list(
  ind_md_ssr       = "Individuals missing data",
  ssr_md_filtering = "Loci missing data",
  null_allele      = "Null alleles",
  LD_filtering     = "Linked loci (LD)"
)
 
filter_params <- list(
  ind_md_ssr       = paste0("MAX_MISSING_IND=", ind_md),
  ssr_md_filtering = paste0("MAX_MISSING_LOCI=", ssr_md),
  null_allele      = paste0("low_quar=", low_quar, ", high_quar=", high_quar, ", IQR_multi=", multi),
  LD_filtering     = paste0("pval_threshold=", threshold_pval)
)
 
# Initial row (before any filtering)
summary_rows[[1]] <- data.frame(
  File   = input_name,
  Step   = 0,
  Filter = "Raw data",
  Params = "-",
  n_ind  = nInd(genind_file),
  n_loci = nLoc(genind_file),
  stringsAsFactors = FALSE
)
 
for (i in seq_along(filters)) {
  step <- filters[i]
 
  if (step == "ind_md_ssr") {
    message("==> Filter: ind_md_ssr (threshold = ", ind_md, ")")
    genind_file <- ind_md_filtering(genind_file, ind_md)
 
  } else if (step == "ssr_md_filtering") {
    message("==> Filter: ssr_md_filtering (threshold = ", ssr_md, ")")
    genind_file <- ssr_md_filtering(genind_file, ssr_md)
 
  } else if (step == "null_allele") {
    message("==> Filter: null_allele (low_quar=", low_quar,
            ", high_quar=", high_quar, ", multi=", multi, ")")
    genind_file <- null_allele(genind_file, low_quar, high_quar, multi)
 
  } else if (step == "LD_filtering") {
    message("==> Filter: LD_filtering (threshold_pval = ", threshold_pval, ")")
    genind_to_genepop(genind_file, "filter_output/genepop_for_LD.txt")
    genind_file <- LD_filtering("filter_output/genepop_for_LD.txt",
                                threshold_pval, genind_file)
  } else {
    warning("Unknown filter: ", step, " — skipped.")
  }
 
  summary_rows[[i + 1]] <- data.frame(
    File   = input_name,
    Step   = i,
    Filter = filter_labels[[step]],
    Params = filter_params[[step]],
    n_ind  = nInd(genind_file),
    n_loci = nLoc(genind_file),
    stringsAsFactors = FALSE
  )
}
 
# Write filtering summary table
summary_table <- do.call(rbind, summary_rows)
output_summary <- paste0("summary/filtering_summary.txt")
write.table(summary_table,
            file      = output_summary,
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)
 
# Write filtered dataset (same format as input)
df_out <- genind2df(genind_file, sep = "/", usepop = FALSE)
df_out <- data.frame(
  Ind = rownames(df_out),
  Pop = pop(genind_file),
  df_out,
  row.names = NULL
)
 
output_file <- paste0("filter_output/", input_name,"_filtered.txt")
write.table(df_out,
            file      = output_file,
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)