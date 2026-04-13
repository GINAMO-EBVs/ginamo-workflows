##### Load packages #####
library(adegenet)   # genind, find.clusters, dapc, scaleGen
library(ade4)       # dudi.pca
library(vcfR)       # read.vcfR, vcfR2genind
library(ggplot2)    # figure PCA
library(dplyr)

##### Load arguments #####
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    input_file    = NULL,
    input_format  = NULL,    # "VCF" | "SSR"
    k_min         = 1L,
    k_max         = 10L,
    k_fixed       = NULL,    # NULL = auto, integer = K fixed
    n_pca_max     = 300L,
    n_start       = 10L,     #only on auto
    seed          = 42L      #NULL = random, integer = reproductible
  )
  
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
           "--input"     = { params$input_file   <- args[i+1]; i <- i+2 },
           "--format"    = { params$input_format  <- args[i+1]; i <- i+2 },
           "--k-min"     = { params$k_min         <- as.integer(args[i+1]); i <- i+2 },
           "--k-max"     = { params$k_max         <- as.integer(args[i+1]); i <- i+2 },
           "--k-fixed"   = { params$k_fixed       <- as.integer(args[i+1]); i <- i+2 },
           "--n-pca-max" = { params$n_pca_max     <- as.integer(args[i+1]); i <- i+2 },
           "--n-start"   = { params$n_start       <- as.integer(args[i+1]); i <- i+2 },
           "--seed"      = {
             val <- args[i+1]
             params$seed <- if (val == "NULL") NULL else as.integer(val)
             i <- i+2
           },
           { stop(paste("Unknown argument :", args[i])) }
    )
  }
  return(params)
}

params <- parse_args(args)

cat("=== dAPC ===\n")
cat("File :", params$input_file, "\n")
cat("Format  :", params$input_format, "\n")
cat("K mode  :", ifelse(is.null(params$k_fixed), "auto", paste("fixed =", params$k_fixed)), "\n\n")

##### Load files + genind conversion #####
cat("[1/5] Loading data and genind conversion...\n")

load_genetic_data <- function(file, format) {
  if (toupper(format) == "VCF") {
    vcf   <- vcfR::read.vcfR(file, verbose = FALSE)
    gdata <- vcfR::vcfR2genind(vcf)
  } else {
    ssr_raw <- read.csv(file, sep="\t")
    geno_cols <- setdiff(colnames(ssr_raw), c("Ind", "Pop"))
    
    gdata <- df2genind(ssr_raw[geno_cols],
                       pop = ssr_raw$Pop,
                       ind.names = ssr_raw$Ind,
                       ploidy = 2,
                       NA.char = c("0", "0/0", "NA","NA/NA"),
                       sep = "/")
  }
  
  cat(sprintf("  -> %d individuals | %d loci | %d allele columns\n",
              nInd(gdata), nLoc(gdata), ncol(gdata@tab)))
  return(gdata)
}

# Load genetic data
gdata <- load_genetic_data(params$input_file, params$input_format)

if (nInd(gdata) < 10) stop("Too few individuals (< 10).")
if (nLoc(gdata) < 5)  stop("Too few loci (< 5).")

##### Data cleaning #####
cat("[1b/5] Cleaning data (monomorphic loci, all-NA columns)...\n")

# 1. Remove loci where ALL individuals are NA
na_prop <- apply(gdata@tab, 2, function(x) mean(is.na(x)))
all_na_cols <- which(na_prop == 1)
if (length(all_na_cols) > 0) {
  # Convert column indices to locus names for dropLoci
  loci_to_drop <- unique(gdata@loc.fac[all_na_cols])
  gdata <- gdata[loc = setdiff(locNames(gdata), as.character(loci_to_drop))]
  cat(sprintf("  -> Removed %d all-NA loci.\n", length(loci_to_drop)))
}

# 2. Remove monomorphic loci (allele freq = 0 or 1 across all individuals,
#    i.e. variance = 0 after NA removal) — these cause Inf/NaN in svd()
tab_nona <- gdata@tab
tab_nona[is.na(tab_nona)] <- 0  # temporary substitution to compute variance
locus_var <- apply(tab_nona, 2, var)
mono_cols <- which(locus_var == 0)
if (length(mono_cols) > 0) {
  loci_mono <- unique(gdata@loc.fac[mono_cols])
  gdata <- gdata[loc = setdiff(locNames(gdata), as.character(loci_mono))]
  cat(sprintf("  -> Removed %d monomorphic loci.\n", length(loci_mono)))
}

cat(sprintf("  -> Retained: %d individuals | %d loci | %d allele columns\n",
            nInd(gdata), nLoc(gdata), ncol(gdata@tab)))

if (nLoc(gdata) < 5) stop("Too few polymorphic loci remaining after cleaning (< 5).")

##### PCA -- scaleGen() + dudi.pca() #####

cat("[2/5] PCA (scaleGen + dudi.pca)...\n")

# Theoretical bound for the rank: min(n-1, p-1) because centring implies a linear constraint
n_pca    <- min(params$n_pca_max, nInd(gdata) - 1L, ncol(gdata@tab) - 1L)
X_scaled <- scaleGen(gdata, NA.method = "mean", scale = FALSE, center = TRUE)

pca_result <- dudi.pca(
  X_scaled,
  cent   = FALSE, 
  scale  = FALSE,
  scannf = FALSE,
  nf     = n_pca
)

# Actual rank returned by dudi.pca (may differ by 1 from the theorical bound according to ade4)
# n_pca is overwritten to ensure consistency across all subsequent steps
n_pca <- length(pca_result$eig)

var_explained <- pca_result$eig / sum(pca_result$eig) * 100
cum_var        <- cumsum(var_explained)

cat(sprintf("  -> %d calculated PC axes (effective rank).\n", n_pca))
cat(sprintf("  -> Cumulative variance across %d axes : %.1f%%\n", n_pca, cum_var[n_pca]))

pca_scores <- pca_result$li  # n x n_pca matrix used for clustering

##### K INFERENCE #####

cat("[3/5] K inference...\n")

n <- nrow(pca_scores)
p <- ncol(pca_scores)

if (!is.null(params$k_fixed)) { # Fixed mode: K set by user -- inference is skipped
  k_infer    <- params$k_fixed
  k_adegenet <- NA_integer_
  k_bic      <- NA_integer_
  cat(sprintf("  -> K fixed by user : %d\n", k_infer))
  
} else { #automatic mode : BIC k-means
  k_range    <- params$k_min:params$k_max
  bic_scores <- rep(NA_real_, length(k_range))
  
  # find.clusters adegenet (BIC)
  cat(" find.clusters adegenet (internal BIC)...\n")
  if (!is.null(params$seed)) set.seed(params$seed)
  fc_auto <- tryCatch(
    find.clusters(
      gdata,
      n.pca       = n_pca,
      n.clust     = NULL,
      method      = "kmeans",
      stat        = "BIC",
      max.n.clust = params$k_max,
      n.start     = params$n_start,
      choose      = FALSE,
      graph       = FALSE
    ),
    error = function(e) { warning(paste("find.clusters failed :", e$message)); NULL }
  )
  
  k_adegenet <- if (!is.null(fc_auto)) {
    as.integer(sub("K=", "", names(which.min(fc_auto$Kstat))))
  } else {
    NA_integer_
  }
  cat(sprintf("      -> K adegenet : %s\n",
              ifelse(is.na(k_adegenet), "failed", k_adegenet)))
  
  k_infer <- k_adegenet
}

# Final k-means clustering at the retained K
if (!is.null(params$seed)) set.seed(params$seed)
km_final <- kmeans(pca_scores, centers = k_infer,
                   iter.max = 1e4, nstart = params$n_start)
fc_final <- list(grp = factor(km_final$cluster))

cat("\n  -> Group sizes:\n")
print(table(fc_final$grp))

##### DAPC -- paxes = k-1 (Thia 2022 rule) #####

# Minimum 1 PCA axis; cap to the effective PCA rank to avoid svd() overflow
if (k_infer > 1) {
  n_paxes <- max(1L, min(k_infer - 1L, n_pca))
  n_daxes <- max(1L, k_infer - 1L)
} else {
  stop(paste0("\n[dAPC Error]: k_infer = ", k_infer, 
              ". Discriminant Analysis requires at least 2 clusters.\n",
              "The data suggests a panmictic population (no genetic structure)."))
}


cat(sprintf("[4/5] DAPC (paxes = %d, n.da = %d)...\n", n_paxes, n_daxes))

dapc_result <- dapc(
  gdata,
  pop   = fc_final$grp,
  n.pca = n_paxes,
  n.da  = n_daxes,
  graph = FALSE
)

assign_success <- mean(dapc_result$assign == fc_final$grp) * 100
cat(sprintf("  -> Reassignment rate : %.1f%%\n", assign_success))

##### Exports files #####

cat("[5/5] Exporting results...\n")

# Scores LD
ld_scores <- cbind(
  Individual = rownames(dapc_result$ind.coord),
  Cluster_kmeans = as.character(fc_final$grp),
  as.data.frame(dapc_result$ind.coord)
)
write.table(ld_scores, "outputs/output_ld.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Assignments + membership probabilities
assign_df <- data.frame(
  Individual     = rownames(dapc_result$ind.coord),
  Cluster_kmeans = as.character(fc_final$grp),
  Cluster_dapc   = as.character(dapc_result$assign),
  Match          = dapc_result$assign == fc_final$grp,
  as.data.frame(dapc_result$posterior)
)
write.table(assign_df, "outputs/output_assign.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary stats
adn_line <- paste("K adegenet [kmeans BIC] :",
                  ifelse(is.na(k_adegenet), "N/A (fixed mode)", k_adegenet))

writeLines(c(
  "=== Automated DAPC Summary ===",
  paste("Date               :", Sys.time()),
  paste("File               :", params$input_file),
  paste("Format             :", params$input_format),
  paste("Seed               :", params$seed),
  paste("N individuals      :", nInd(gdata)),
  paste("N loci             :", nLoc(gdata)),
  paste("N total alleles    :", ncol(gdata@tab)),
  paste("N PCA axes         :", n_pca),
  paste("Cumulative var.    :", sprintf("%.1f%%", cum_var[n_pca])),
  "",
  "=== K Inference ===",
  adn_line,
  paste("K retained         :", k_infer),
  paste("Mode               :", ifelse(is.null(params$k_fixed), "auto", "fixed")),
  "",
  "=== DAPC ===",
  paste("paxes (k-1)        :", n_paxes),
  paste("DA axes            :", n_daxes),
  paste("Reassignment       :", sprintf("%.2f%%", assign_success)),
  "",
  "=== Group sizes ===",
  capture.output(print(table(fc_final$grp)))
), "outputs/output_stats.txt")

##### Shared colour palette (used across all figures) #####
pal_anchors  <- c("#7B2D8B", "#2166AC", "#41B6C4", "#1A9641",
                  "#D9EF8B", "#FD8D3C", "#D73027")
cluster_cols <- colorRampPalette(pal_anchors)(k_infer)

##### BIC curve (auto mode only) #####
if (is.null(params$k_fixed) && !is.null(fc_auto)) {
  bic_df <- data.frame(
    K   = as.integer(sub("K=", "", names(fc_auto$Kstat))),
    BIC = as.numeric(fc_auto$Kstat)
  )
  
  bic_plot <- ggplot(bic_df, aes(x = K, y = BIC)) +
    geom_line(colour = "grey50", linewidth = 0.7) +
    geom_point(size = 3, colour = "grey30") +
    geom_point(data = subset(bic_df, K == k_infer),
               size = 4, colour = "#D73027", shape = 18) +
    geom_vline(xintercept = k_infer, linetype = "dashed",
               colour = "#D73027", linewidth = 0.5) +
    scale_x_continuous(breaks = bic_df$K) +
    labs(
      title    = "K-means BIC — K selection",
      subtitle = paste0("Retained K = ", k_infer, "  (red diamond)"),
      x = "Number of clusters K",
      y = "BIC"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("outputs/output_bic.png", plot = bic_plot,
         width = 7, height = 5, dpi = 150, bg = "white")
  cat("  -> BIC curve     :", "outputs/output_bic.png", "\n")
} else {
  cat("  -> BIC curve     : skipped (fixed K mode)\n")
}

##### Scree plot (PCA variance explained) #####
# Show at most the first 30 axes for readability
n_show     <- min(30L, n_pca)
scree_df   <- data.frame(
  Axis      = seq_len(n_show),
  Var       = var_explained[seq_len(n_show)],
  CumVar    = cum_var[seq_len(n_show)]
)

scree_plot <- ggplot(scree_df, aes(x = Axis)) +
  geom_col(aes(y = Var), fill = "#2166AC", alpha = 0.75, width = 0.7) +
  geom_line(aes(y = CumVar), colour = "#D73027",
            linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = CumVar), colour = "#D73027", size = 1.8) +
  scale_x_continuous(breaks = seq_len(n_show)) +
  labs(
    title    = "PCA scree plot",
    subtitle = sprintf("First %d axes shown  |  total rank = %d", n_show, n_pca),
    x = "PC axis",
    y = "Variance explained (%)"
  ) +
  # Secondary axis label for cumulative variance (cosmetic only)
  annotate("text", x = n_show * 0.92, y = max(scree_df$CumVar) * 0.97,
           label = "— cumulative", colour = "#D73027", size = 3.2, hjust = 1) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 7))

ggsave("outputs/output_scree.png", plot = scree_plot,
       width = 8, height = 5, dpi = 150, bg = "white")
cat("  -> Scree plot    :", "outputs/output_scree.png", "\n")

##### PCA colored by K-means cluster #####
pca_df <- data.frame(
  PC1     = pca_result$li[, 1],
  PC2     = pca_result$li[, 2],
  Cluster = factor(as.character(fc_final$grp),
                   levels = sort(unique(as.character(fc_final$grp))))
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Cluster, fill = Cluster)) +
  stat_ellipse(geom = "polygon", level = 0.67, alpha = 0.12, linewidth = 0.4) +
  geom_point(size = 2.2, alpha = 0.85, shape = 21,
             colour = "white", stroke = 0.3, aes(fill = Cluster)) +
  scale_colour_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  labs(
    title    = paste0("PCA - K-means clusters (K = ", k_infer, ")"),
    subtitle = ifelse(is.null(params$k_fixed),
                      paste0("K adegenet = ", ifelse(is.na(k_adegenet), "?", k_adegenet),
                             "  |  retained = ", k_infer),
                      paste0("K fixed by user = ", k_infer)),
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

ggsave("outputs/output_pca.png", plot = pca_plot,
       width = 8, height = 6, dpi = 150, bg = "white")
cat("  -> PCA figure:", "outputs/output_pca.png", "\n")

##### DAPC scatter (DA axes) #####
dapc_coords <- as.data.frame(dapc_result$ind.coord)
dapc_coords$Cluster <- factor(as.character(fc_final$grp),
                              levels = sort(unique(as.character(fc_final$grp))))

if (ncol(dapc_result$ind.coord) >= 2) {
  # Two or more DA axes: standard 2-D scatter
  dapc_scatter <- ggplot(dapc_coords,
                         aes(x = LD1, y = LD2, colour = Cluster, fill = Cluster)) +
    stat_ellipse(geom = "polygon", level = 0.67, alpha = 0.12, linewidth = 0.4) +
    geom_point(size = 2.2, alpha = 0.85, shape = 21,
               colour = "white", stroke = 0.3, aes(fill = Cluster)) +
    scale_colour_manual(values = cluster_cols) +
    scale_fill_manual(values = cluster_cols) +
    labs(
      title    = paste0("DAPC scatter (K = ", k_infer, ")"),
      subtitle = sprintf("Reassignment rate: %.1f%%", assign_success),
      x = "DA axis 1",
      y = "DA axis 2"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey40"))
} else {
  # Only one DA axis (K = 2): density plot per cluster
  dapc_scatter <- ggplot(dapc_coords, aes(x = LD1, fill = Cluster, colour = Cluster)) +
    geom_density(alpha = 0.4, linewidth = 0.6) +
    scale_colour_manual(values = cluster_cols) +
    scale_fill_manual(values = cluster_cols) +
    labs(
      title    = paste0("DAPC density — DA axis 1 (K = ", k_infer, ")"),
      subtitle = sprintf("Reassignment rate: %.1f%%", assign_success),
      x = "DA axis 1",
      y = "Density"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey40"))
}

ggsave("outputs/output_dapc_scatter.png", plot = dapc_scatter,
       width = 8, height = 6, dpi = 150, bg = "white")
cat("  -> DAPC scatter  :", "outputs/output_dapc_scatter.png", "\n")

##### Posterior membership barplot (STRUCTURE-style) #####
post_df <- as.data.frame(dapc_result$posterior)
post_df$Individual <- rownames(post_df)
post_df$Cluster    <- factor(as.character(fc_final$grp),
                             levels = sort(unique(as.character(fc_final$grp))))

# Sort individuals by cluster then by their dominant posterior probability
post_df <- post_df[order(post_df$Cluster,
                         -apply(dapc_result$posterior, 1, max)), ]
post_df$Individual <- factor(post_df$Individual, levels = post_df$Individual)

# Reshape to long format without reshape2/tidyr
post_long <- do.call(rbind, lapply(seq_len(k_infer), function(k) {
  col_name <- colnames(dapc_result$posterior)[k]
  data.frame(
    Individual  = post_df$Individual,
    Cluster_ref = post_df$Cluster,
    Component   = col_name,
    Probability = post_df[[col_name]]
  )
}))
post_long$Component <- factor(post_long$Component,
                              levels = colnames(dapc_result$posterior))

barplot_post <- ggplot(post_long,
                       aes(x = Individual, y = Probability, fill = Component)) +
  geom_col(width = 1, colour = NA) +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  # Vertical separators between clusters
  geom_vline(
    xintercept = cumsum(as.numeric(table(post_df$Cluster))) + 0.5,
    colour = "white", linewidth = 0.6
  ) +
  labs(
    title    = paste0("Posterior membership probabilities (K = ", k_infer, ")"),
    subtitle = sprintf("Individuals sorted by cluster then by dominant probability  |  reassignment: %.1f%%",
                       assign_success),
    x = NULL,
    y = "Membership probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank()
  )

# Scale width with number of individuals (min 8, max 24 inches)
fig_width <- max(8, min(24, nInd(gdata) / 15))
ggsave("outputs/output_posterior.png", plot = barplot_post,
       width = fig_width, height = 4, dpi = 150, bg = "white")
cat("  -> Posterior bar :", "outputs/output_posterior.png", "\n")

cat("\n=== Analysis completed successfully ===\n")

