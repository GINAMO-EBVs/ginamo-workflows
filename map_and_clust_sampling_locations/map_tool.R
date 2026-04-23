##### Load packages #####
library(dplyr)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggforce)
library(dbscan)


##### Load arguments #####
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    input_file        = NULL, 
    lon_col           = "Longitude",
    lat_col           = "Latitude",
    cluster_mode      = "auto",     # "auto" | "manual"
    cluster_col       = NULL,
    buffer_km         = NULL,         # buffer radius in km (automatic mode)
    view_xmin         = NULL,        # manual bbox: longitude min
    view_xmax         = NULL,        # manual bbox: longitude max
    view_ymin         = NULL,        # manual bbox: latitude min
    view_ymax         = NULL,        # manual bboxlle : latitude max
    view_countries    = NULL,        # countries to display
    map_title         = NULL,
    width             = 12,
    height            = 9,
    dpi               = 300L
  )
  
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
           "--input"            = { params$input_file       <- args[i+1];              i <- i+2 },
           "--lon-col"          = { params$lon_col          <- args[i+1];              i <- i+2 },
           "--lat-col"          = { params$lat_col          <- args[i+1];              i <- i+2 },
           "--cluster-mode"     = { params$cluster_mode     <- args[i+1];              i <- i+2 },
           "--cluster-col"      = { params$cluster_col      <- args[i+1];              i <- i+2 },
           "--buffer-km"        = { params$buffer_km        <- as.numeric(args[i+1]); i <- i+2 },
           "--show-hull"        = { params$show_hull        <- as.logical(args[i+1]); i <- i+2 },
           "--view-xmin"        = { params$view_xmin        <- as.numeric(args[i+1]); i <- i+2 },
           "--view-xmax"        = { params$view_xmax        <- as.numeric(args[i+1]); i <- i+2 },
           "--view-ymin"        = { params$view_ymin        <- as.numeric(args[i+1]); i <- i+2 },
           "--view-ymax"        = { params$view_ymax        <- as.numeric(args[i+1]); i <- i+2 },
           "--view-countries"   = { params$view_countries   <- args[i+1];              i <- i+2 },
           "--map-title"        = { params$map_title        <- args[i+1];              i <- i+2 },
           { stop(paste("Unknown argument :", args[i])) }
    )
  }
  return(params)
}

params <- parse_args(args)

params$map_padding       = 0.1       # margin around points (extended fraction) — ignored if view_* is provided

##### Check inputs #####

if (!file.exists(params$input_file)) {
  stop(paste0("ERROR: File not found:  ", params$input_file))
}

if (params$cluster_mode == "manual" && is.null(params$cluster_col)) {
  stop("ERROR: Cluster column name is required in 'manual' mode.")
}

##### Read data #####
df <- tryCatch(
  {
    data <- read.table(params$input_file,
                       header    = TRUE,
                       sep       = "\t",
                       quote     = "",
                       fill      = TRUE,
                       comment.char = "",
                       stringsAsFactors = FALSE,
                       check.names = FALSE
    )
    
    colnames(data) <- gsub('"', '', colnames(data))
    
    data
  },
  error = function(e) stop("ERROR while reading the file : ", e$message)
)

# Check mandatory columns
for (col in c(params$lon_col, params$lat_col)) {
  if (!col %in% colnames(df)) {
    stop(paste0("ERROR: Column not found in the file : '", col, "'"))
  }
}

# Rename columns longitude and latitude
df <- df %>% rename(
  .lon = !!sym(params$lon_col),
  .lat = !!sym(params$lat_col)
)

# Conversion
df$.lon <- as.numeric(df$.lon)
df$.lat <- as.numeric(df$.lat)

# Remove missing coords
n_before <- nrow(df)
df <- df %>% filter(!is.na(.lon) & !is.na(.lat))
n_removed <- n_before - nrow(df)
if (n_removed > 0) {
  message("  >> ", n_removed, " line(s) ignored (missing coordinates).")
}
if (nrow(df) == 0) {
  stop("ERROR: No valid rows after filtering for missing coordinates.")
}

##### Clustering #####

message(">> Clustering (mode : ", params$cluster_mode, ")...")

if (params$cluster_mode == "manual") {
  if (!params$cluster_col %in% colnames(df)) {
    stop(paste0("ERROR : column cluster '", params$cluster_col, "' not found."))
  }
  df$.cluster <- as.character(df[[params$cluster_col]])
  message("  >> ", length(unique(df$.cluster)), " cluster(s) detected from the column'", params$cluster_col, "'.")
  
} else {
  # Mode auto : buffer overlap
  # Two individuals belong to the same population if their buffers overlap,
  # i.e. if the point-to-point distance is less than or equal to 2 * buffer_km.
  
  threshold_m <- params$buffer_km * 2 * 1000
  
  sf_pts   <- st_as_sf(df, coords = c(".lon", ".lat"), crs = 4326)
  sf_pts_m <- st_transform(sf_pts, crs = 3857)
  coords   <- st_coordinates(sf_pts_m)
  
  db_res <- dbscan::dbscan(coords, eps = threshold_m, minPts = 1)
  
  # Attribution des IDs (dbscan commence à 1, pas de 0/bruit ici car minPts=1)
  df$.cluster <- paste0("Pop_", db_res$cluster)
}

##### Map window #####
# Priority : manual bbox > countries list > auto (points + padding)

bbox_mode <- "auto"

if (!is.null(params$view_xmin) && !is.null(params$view_xmax) &&
    !is.null(params$view_ymin) && !is.null(params$view_ymax)) {
  
  # Mode 1 : manual bbox provided by the user
  bbox_mode <- "manual"
  bbox <- c(
    left   = params$view_xmin,
    bottom = params$view_ymin,
    right  = params$view_xmax,
    top    = params$view_ymax
  )
  message(">> Map window (manual mode) : lon [", bbox["left"], ", ", bbox["right"],
          "] / lat [", bbox["bottom"], ", ", bbox["top"], "]")
  
} else if (!is.null(params$view_countries)) {
  
  # Mode 2 : countries list — we calculate the union bounding box of their geometries
  bbox_mode <- "countries"
  country_list <- trimws(strsplit(params$view_countries, ",")[[1]])
  message(">> Map window (country view): ", paste(country_list, collapse = ", "))
  
  if (!exists("world")) {
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  }
  
  # Search for countries by name (name, name_long, admin, iso_a2, iso_a3)
  matched <- world[
    tolower(world$name)      %in% tolower(country_list) |
      tolower(world$name_long) %in% tolower(country_list) |
      tolower(world$admin)     %in% tolower(country_list) |
      toupper(world$iso_a2)    %in% toupper(country_list) |
      toupper(world$iso_a3)    %in% toupper(country_list),
  ]
  
  not_found <- country_list[
    !tolower(country_list) %in% tolower(c(world$name, world$name_long, world$admin)) &
      !toupper(country_list) %in% toupper(c(world$iso_a2, world$iso_a3))
  ]
  if (length(not_found) > 0) {
    warning("Country(ies) not found in the Natural Earth database:",
            paste(not_found, collapse = ", "),
            ". Check the names (in English) or ISO codes.")
  }
  
  if (nrow(matched) == 0) {
    warning("No countries recognised, returning to automatic mode.")
    bbox_mode <- "auto"
  } else {
    # Some countries (France, Spain, Portugal, the Netherlands, etc.) 
    #have overseas territories included in their Natural Earth geometry, 
    #which causes the bounding box to blow up
    # Solution: for each entity, the MULTIPOLYGON is broken down into individual 
    # polygons, and only those whose centroid lies within the bounding box of 
    # the largest polygon (mainland territory) are retained
    filter_continental <- function(sf_row) {
      geom <- sf::st_geometry(sf_row)[[1]]
      
      #Break down into individual polygons
      parts <- sf::st_cast(sf::st_sfc(geom, crs = 4326), "POLYGON")
      if (length(parts) == 1) return(sf_row)  # déjà un seul polygone
      
      # Keep only the bigger polygone
      areas  <- as.numeric(sf::st_area(parts))
      main   <- parts[which.max(areas)]
      bb_main <- sf::st_bbox(main)
      
      # Filter: keep the parts whose centroid lies within ±30° of the main one
      centroids_x <- sapply(seq_along(parts), function(i)
        sf::st_coordinates(sf::st_centroid(parts[i]))[1])
      centroids_y <- sapply(seq_along(parts), function(i)
        sf::st_coordinates(sf::st_centroid(parts[i]))[2])
      
      keep <- which(
        centroids_x >= (bb_main["xmin"] - 30) &
          centroids_x <= (bb_main["xmax"] + 30) &
          centroids_y >= (bb_main["ymin"] - 30) &
          centroids_y <= (bb_main["ymax"] + 30)
      )
      filtered_geom <- sf::st_union(parts[keep])
      sf::st_geometry(sf_row) <- sf::st_sfc(filtered_geom, crs = 4326)
      sf_row
    }
    
    matched_cont <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i)
      filter_continental(matched[i, ])
    ))
    
    bb <- sf::st_bbox(matched_cont)
    
    # Merge with the user's point bbox :
    # ensures that all points remain visible even if some are outside the 
    # mainland (e.g. points in French Guiana where `view_countries` is set to "France")..
    bbox <- c(
      left   = min(as.numeric(bb["xmin"]) - 0.5, min(df$.lon) - 0.5),
      bottom = min(as.numeric(bb["ymin"]) - 0.5, min(df$.lat) - 0.5),
      right  = max(as.numeric(bb["xmax"]) + 0.5, max(df$.lon) + 0.5),
      top    = max(as.numeric(bb["ymax"]) + 0.5, max(df$.lat) + 0.5)
    )
    
    # Notify if any points are outside the continental territory of the requested countries
    pts_out <- df[
      df$.lon < as.numeric(bb["xmin"]) | df$.lon > as.numeric(bb["xmax"]) |
        df$.lat < as.numeric(bb["ymin"]) | df$.lat > as.numeric(bb["ymax"]),
    ]
    if (nrow(pts_out) > 0) {
      message("  >> WARNING : ", nrow(pts_out), " point(s) outside the mainland ",
              "of selected countries — the scope has been expanded to include them.")
    }
    
    message("  >> ", nrow(matched), " countries found. Windows : lon [",
            round(bbox["left"], 2), ", ", round(bbox["right"], 2),
            "] / lat [", round(bbox["bottom"], 2), ", ", round(bbox["top"], 2), "]")
  }
}


if (bbox_mode == "auto") {
  # Mode 3 : automatic —  points + padding
  lon_range <- range(df$.lon)
  lat_range <- range(df$.lat)
  lon_pad   <- max((lon_range[2] - lon_range[1]) * params$map_padding, 0.5)
  lat_pad   <- max((lat_range[2] - lat_range[1]) * params$map_padding, 0.5)
  bbox <- c(
    left   = lon_range[1] - lon_pad,
    bottom = lat_range[1] - lat_pad,
    right  = lon_range[2] + lon_pad,
    top    = lat_range[2] + lat_pad
  )
  message(">> Map window (automatic mode) : lon [", round(bbox["left"], 4), ", ",
          round(bbox["right"], 4), "] / lat [", round(bbox["bottom"], 4), ", ",
          round(bbox["top"], 4), "]")
}

##### Vector base map (rnaturalearth) #####
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
message("  >> Vector base map load (", nrow(world), " countries).")

##### Colour palette #####
clusters      <- sort(unique(df$.cluster))
n_clusters    <- length(clusters)
if (n_clusters <= 12) {
  palette <- RColorBrewer::brewer.pal(max(3, n_clusters), "Set1")[seq_len(n_clusters)]
} else {
  palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
}
names(palette) <- clusters

##### Build map #####
message(">> Building map...")


# Base : country outlines (rnaturalearth)
p <- ggplot() +
  geom_sf(
    data  = world,
    fill  = "grey92",
    colour = "grey60",
    linewidth = 0.3
  ) +
  theme_bw(base_size = 12)


# Points
p <- p +
  geom_point(
    data    = df,
    mapping = aes(x = .lon, y = .lat, colour = .cluster), 
    shape = 16,
    size    = 2.5,
    stroke  = 0.7
  ) +
  scale_colour_manual(values = palette, name = "Population") +
  scale_fill_manual(values   = palette, name = "Population")

# Ellipses by population (added after the map was finalised)
  # Calculation of ellipse parameters by cluster (centre + semi-axes + angle)
ellipse_params <- do.call(rbind, lapply(clusters, function(cl) {
  pts <- df[df$.cluster == cl, c(".lon", ".lat")]
  n   <- nrow(pts)
  if (n == 1) {
    # 1 point : small fixed circle
    data.frame(
      .cluster = cl,
      x0 = pts$.lon[1], y0 = pts$.lat[1],
      a  = 0.3, b = 0.3, angle = 0
    )
    
  } else if (n == 2) {
    # 2 points : an ellipse elongated along the line segment
    cx    <- mean(pts$.lon)
    cy    <- mean(pts$.lat)
    dx    <- diff(pts$.lon)
    dy    <- diff(pts$.lat)
    angle <- atan2(dy, dx)
    
    a     <- sqrt(dx^2 + dy^2) / 2 + 0.3   # semi-major axis
    b     <- 0.3                              # fixed minor semi-axis
  
    data.frame(.cluster = cl, x0 = cx, y0 = cy, a = a, b = b, angle = angle)

  } else {
    # 3+ points: ellipse based on the standard deviation (± 2 SD = ~95% of the points)
    cx    <- mean(pts$.lon)
    cy    <- mean(pts$.lat)
    # ACP to orient the ellipse along the main direction of the cloud
    pca   <- prcomp(pts, center = TRUE, scale. = FALSE)
    angle <- atan2(pca$rotation[2, 1], pca$rotation[1, 1])
    # Standard deviations in the PCA's own coordinate system
    scores <- pca$x
    a     <- max(sd(scores[, 1]) * 2, 0.15)
    b     <- max(sd(scores[, 2]) * 2, 0.15)
    data.frame(.cluster = cl, x0 = cx, y0 = cy, a = a, b = b, angle = angle)
  }
}))

  # Semi-transparent filling
  p <- p +
    ggforce::geom_ellipse(
      data    = ellipse_params,
      mapping = aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle,
                    fill = .cluster),
      alpha   = 0.18,
      colour  = NA
    ) +
    ggforce::geom_ellipse(
      data    = ellipse_params,
      mapping = aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle,
                    colour = .cluster),
      fill     = NA,
      linewidth = 0.55,
      linetype  = "dashed"
  )

# Card design
p <- p +
  labs(
    title    = params$map_title,
    subtitle = paste0(nrow(df), " samples — ", n_clusters, " population(s) — Buffer : ",
                      if (params$cluster_mode == "auto") paste0(params$buffer_km, " km") else "manual mode"),
    x        = "Longitude",
    y        = "Latitude",
    caption  = paste0(
      "Geographical reference system : WGS84",
      " | Base map : Natural Earth",
      " | Generated : ", format(Sys.Date(), "%Y-%m-%d")
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    plot.caption     = element_text(size = 7,  colour = "grey55", hjust = 1),
    axis.text        = element_text(size = 9),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3)
  ) +
  # coord_sf declared last to force the window to apply to all geom_sf
  coord_sf(
    xlim   = c(bbox["left"],   bbox["right"]),
    ylim   = c(bbox["bottom"], bbox["top"]),
    expand = FALSE
  )

#Save map
ggsave("outputs/samples_map.png", plot = p, width = params$width, height = params$height,
         dpi = params$dpi, device = "png")
message("  >> Saved map.")

# Export table with new pops

# Reconstruction of the output table
if (params$cluster_mode=="auto"){
  df_out <- df %>%
    rename(
      !!params$lon_col := .lon,
      !!params$lat_col := .lat,
      Population    := .cluster
    ) %>%
    select(-any_of(".id"))
  
  
  write.table(df_out, file = "outputs/table_output.txt", sep = "\t",
              quote = FALSE, row.names = FALSE)
}

