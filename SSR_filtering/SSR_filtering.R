# Load packages
library(adegenet)
library(poppr)
library(hierfstat)
library(genepop)

#Load arguments
args <- commandArgs(trailingOnly = TRUE)

input_path    <- args[1]   # path to input file
input_name <- basename(args[2])
# Extract base name (remove extension, handle parentheses from Galaxy labels)
if (grepl("\\([^)]+\\)\\s*$", input_name)) {
  input_name <- sub(".*\\(([^)]+)\\)\\s*$", "\\1", input_name)
} else {
  input_name <- sub("\\.[^.]+$", "", input_name)
}  

ind_md <- as.numeric(args[3])   # NA if not selected
order_ind <- as.numeric(args[4]) # NA if not selected
ssr_md <- as.numeric(args[5])   # NA if not selected
order_ssr <- as.numeric(args[6]) # NA if not selected
low_quar <- as.numeric(args[7])   # NA if not selected
high_quar <- as.numeric(args[8])   # NA if not selected
multi <- as.numeric(args[9])   # NA if not selected
order_null <- as.numeric(args[10]) # NA if not selected
threshold_pval <- as.numeric(args[11])  # NA if not selected
order_LD <- as.numeric(args[12]) # NA if not selected

# Create a list with the order to applied the different filters
filter_list <- list(
  list(name = "ind_md_ssr",       order = order_ind,  active = !is.na(order_ind)),
  list(name = "ssr_md_filtering", order = order_ssr,  active = !is.na(order_ssr)),
  list(name = "null_allele",      order = order_null, active = !is.na(order_null)),
  list(name = "LD_filtering",     order = order_ld,   active = !is.na(order_ld))
)

active_filters <- Filter(function(f) f$active, filter_list)
active_filters <- active_filters[order(sapply(active_filters, function(f) f$order))]
filters <- sapply(active_filters, function(f) f$name)

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
  
  # Détection automatique du digit_num via alleles()
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


########################################################################################
# Function : LD_filtering
# Description : Performs the exact test for each pair of loci for each population. 
#               Deletes the locus involved  in the most pairs with significant p-values.
#########################################################################################

LD_filtering <- function(genpop_file_path, threshold_pval, genind_file) { 
  
  #Exact test from genepop package
  test_LD(genepop_file_path, outputFile = "LD_test.txt")

  # Read output file from test_LD
  lines <- readLines("LD_test.txt")
  lines <- trimws(lines, which = "right")

  # Extraction of the pairs with significant p-values
    # Extract only the lines containing the results for locus pairs
  global_lines <- lines[grepl("&", lines)]

  global_table <- do.call(rbind, lapply(global_lines, function(l) {
    parts <- strsplit(trimws(l), "\\s+")[[1]]
    #Extract p-values from the 6th position and remove special characters
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

  # Greedy strategy : Deletes the locus involved  in the most pairs with significant p-values.
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
  genind_nLD <- genind_file[loc = loci_a_garder]

  return(genind_nLD)

}

pops <- pop(genind_nLD)
  df <- genind2df(genind_nLD, sep = "/", usepop = FALSE)

  df <- data.frame(
    Ind = rownames(df),
    Pop = pop(genind_nLD),
    df,
    row.names = NULL
  )

  return(df)


########################
# Main execution
# Conversion et test LD sur le genind final (après filtrage missing data et allèles nuls)
genind_to_genepop(genind_sans_outliers, "genepop_SSR_filtering.txt")
##########################

# Summary table
summary_rows <- list()

for (step in filters) {

  nSSR_before <- nLoc(genind_file)
  nInd_before <- nInd(genind_file)

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

  nSSR_after <- nLoc(genind_file)
  nInd_after <- nInd(genind_file)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Step       = step,
    nSSRbefore = nSSR_before,
    nSSRafter  = nSSR_after,
    nIndbefore = nInd_before,
    nIndafter  = nInd_after,
    stringsAsFactors = FALSE
  )
}

# Write filtering summary table
summary_table <- do.call(rbind, summary_rows)
write.table(summary_table,
            file      = "summary/filtering_summary.txt",
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
