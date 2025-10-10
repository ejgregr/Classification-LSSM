#----------------------------------------------------------------------------
# Script:  Broughton2_functions.R
# Created: May 2024 from the original Broughton script. EJG
#
# Purpose: Support building and evaluation of Seaweed Clusters for Coastal BC.
# Illustrated at 3 spatial extents using the Broughton Region.
#
# Notes:
#  - 2024/06/05: Updated from the almost final DFO version.

#================================== Load require packages =================================

# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "reshape2", "tidyr","dplyr", "raster", "stringr", "rasterVis",
                        "RColorBrewer", "factoextra", "ggpubr", "cluster", "rmarkdown","lubridate",
                        "knitr", "tinytex", "kableExtra", "e1071",
                        "ncdf4", "terra", "akima", "sf",
                        "clue", "patchwork" )

# viridisLite an option for gradient palettes

# "diffeR", "vegan", "ranger", "e1071", "forcats", "measures", "caret", "PresenceAbsence"
# "randomForest", "spatialEco", "xlsx", "robustbase", "biomod2", "sp", "magrittr", "tinytex", "rmarkdown", "binr", 'gwxtab'

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)


#===================================== Functions =========================================
# ---- CLUSTER matching function: Align kmean cluster labels by membership overlap 
# Will use Hungarian algorithm. Needs packages("clue"). Workflow by ChatGPT

# Pass the two kmeans to align: km_ref = reference, km_new = to be aligned
# Returns: $mapping (index = new label, value = ref label),
#          $aligned (new labels remapped to ref),
#          $contingency (ref x new table)
AlignClusters <- function(km_ref, km_new) {
  stopifnot(length(km_ref$cluster) == length(km_new$cluster))
  K_ref <- nrow(km_ref$centers); K_new <- nrow(km_new$centers)
  stopifnot(K_ref == K_new)  # same K
  K <- K_ref
  
  ref <- as.integer(km_ref$cluster)
  new <- as.integer(km_new$cluster)
  
  # Contingency: rows = ref labels, cols = new labels
  M <- table(factor(ref, levels = 1:K),
             factor(new, levels = 1:K))
  M <- as.matrix(M)
  
  # Hungarian solves a min-cost; we want to maximize overlap
  cost <- max(M) - M
  assign <- clue::solve_LSAP(cost)  # for each ref row i, chosen new col j
  
  # Build mapping new->ref (index = new label, value = ref label)
  mapping <- integer(K)
  for (i in 1:K) {
    j <- assign[i]
    mapping[j] <- i
  }
  
  aligned <- mapping[new]
  list(mapping = mapping, aligned = aligned, contingency = M)
}

# mapping: named character or integer vector like c("1"="2","2"="3","3"="1")
# Pass: the kmeans to modify and the mapping (as a perm vector from AlignClusters)
# Returns: updated kmeans object
ApplyAlignment <- function(km, perm) {
  stopifnot(inherits(km, "kmeans"))
  k <- nrow(km$centers)
  stopifnot(length(perm) == k, setequal(perm, seq_len(k)))  # must be a permutation
  
  inv <- integer(k); inv[perm] <- seq_len(k)  # inverse: inv[new] = old
  
  km2 <- km
  km2$cluster  <- perm[km$cluster]           # relabel assignments
  km2$centers  <- km$centers[inv, , drop = FALSE]  # reorder to 1..k in new label order
  km2$withinss <- km$withinss[inv]
  km2$size     <- km$size[inv]
  rownames(km2$centers) <- as.character(seq_len(k))
  km2
}

# Given a correlation matrix, return the highest correlation and the variable pair.
# NOTE: Algorithm fails down at nrom(idx) when a single correlation is found in the cortable.
#     ---> Soln is just to lower the threshold. :\
showHighestPairs <- function(preds, threshold) {
  
  cor_table <- round( cor( preds ), 3)
  cor_table[lower.tri(cor_table, diag=TRUE)] <- NA
  
  high_rows <- apply(cor_table, 1, function(row) any(abs(row) > threshold, na.rm = TRUE))
  ctab <- cor_table[ high_rows, ]
  
  ctab2 <- abs(ctab)
  idx <- which(ctab2 >= threshold, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    message("No correlations above threshold.")
    return(invisible(NULL))
  }
  
  for (i in seq_len(nrow(idx))) {
    row_name <- rownames(ctab)[idx[i, 1]]
    col_name <- colnames(ctab)[idx[i, 2]]
    value <- ctab[idx[i, 1], idx[i, 2]]
    cat(sprintf("Row: %s, Col: %s, Correlation: %.3f\n", row_name, col_name, value))
  }
}

#---- Returns a stack of integerized rasters from a raster stack ----
Integerize <- function( in_layers, sig = 1000 ) {
  int_layers <- brick()
  for (i in 1:nlayers(in_layers)) {
    raster_to_scale <- in_layers[[i]]
    int_raster <- as.integer( raster_to_scale * sig )
    int_layers <- stack( int_layers, int_raster )
  }
  
  names( int_layers ) <- names( in_layers )
  return( int_layers )
}

# Uses global variables. Relies on stack_data and transformed stack data (t_stack)
MakeSkewTable <- function() {
  # Target data frame for the skew comparison table
  skdat <- data.frame(Predictor = character(), Pre_Skew = numeric(), Post_Skew = numeric(), 
                      stringsAsFactors = FALSE)
  
  for (i in 1:ncol(stack_data)) {
    # Calculate skewness for each column in both datasets
    pre  <- skewness(stack_data[, i], na.rm = TRUE) 
    post <- skewness(t_stack[, i], na.rm = TRUE)
    
    # Combine the column name and skewness values into a new data frame row
    new_row <- data.frame(
      Predictor = colnames(stack_data)[i],  # Get column name
      Pre_Skew  = round(pre, 3),            # Skewness before transformation
      Post_Skew = round(post, 3)            # Skewness after transformation
    )
    # Append the new row to skdat
    skdat <- rbind(skdat, new_row)
  }
  return(skdat)
  
  # Can hard-code the addition of ceilings and transforms directly to the table 
  # df <- data.frame(
  #   Predictor = c("REI", "freshwater_index", "standard_dev_slope", "temp_range"),
  #   Ceiling = c(0.3, 0.025, 10, NA),
  #   Power_transform = c("1/2", "1/3", "1/2", "1/2")
  # )
  # Merge the data frames
}

# Selected FVCOM data with skews >1:
# salt_DEC_surf_max   CS_DEC_bott_ave   CS_JUL_bott_ave   CS_JUL_surf_ave   CS_JUL_surf_max 
# 1.338474          1.142058          2.048373          1.321944          1.683211 
MakeMoreNormal <- function( dat ){
#  asin(x/max(x))
#  log(x)
  x <- dat[, "salt_DEC_surf_max"]
  y <- x^-3
  dat[, "salt_DEC_surf_max"] <- y
  
  x <- dat[, "CS_DEC_bott_ave"]
  y <- log(x+1)
  dat[, "CS_DEC_bott_ave"] <- y
  
  x <- dat[, "CS_JUL_bott_ave"]
  y <- (x+.1)^-2
  dat[, "CS_JUL_bott_ave"] <- y
  
  x <- dat[, "CS_JUL_surf_ave"]
  y <- (x+.1)^-2
  dat[, "CS_JUL_surf_ave"] <- y
  
  x <- dat[, "CS_JUL_surf_max"]
  y <- (x/max(x)+.5)^-2
  dat[, "CS_JUL_surf_max"] <- y

  #bring lower values in and transform
  x <- dat[, "salt_JUL_surf_min"]
  floor <- 20
  y <- ifelse(x < floor, floor, x)
  y <- (x/max(x)+.5)^-2
  y <- y^-3
  dat[, "salt_JUL_surf_min"] <- y
  
  return( dat ) 
}

PlotHistos <- function( adf, ptitle ){
  df_long <- pivot_longer(adf, everything(), names_to = "variable", values_to = "value")
  ggplot(df_long, aes(x = value)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    facet_wrap(~ variable, scales = "free") +
    theme_minimal() +
    labs(title = ptitle) +
    theme(plot.title = element_text(hjust = 0),
          axis.title.x = element_blank())
}

PlotHeatMap <- function( clust, dat, main_txt = "" ){
  # Define color palette
  pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette
  
  profile_data <- data.frame(
    cluster = factor(clust, levels = seq_len(nclust)),
    dat,
    check.names = FALSE
  )
  
  cluster_sd <- profile_data %>%
    group_by(cluster) %>%
    summarise_all(sd)
  
  x <- as.data.frame( cluster_sd )
  xm <- melt( x, id.var = "cluster" )
  
  z_heat <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
    geom_tile() +
    scale_fill_gradientn(colours = pal_heat) +
    theme_minimal() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    labs(title = main_txt, x = "Clusters", y = "Predictors", fill = "Value")
  return( z_heat )
}

PlotViolins <- function( clusts, dat, pal ){
  x <- dat
  x$cluster <- as.factor( clusts )
  y <- x %>%
    pivot_longer(cols = -cluster, names_to = "predictor", values_to = "value")
  
  # Create violin plot
  vplots <- 
    ggplot(y, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(trim = FALSE) +
    facet_wrap(~ predictor, scales = "free_y") +
    scale_fill_manual(values = pal) +
    theme_minimal() +
    labs(title = "",
         x = "Cluster",
         y = "Value")
  
  return( vplots )
}

PlotMap <- function(coords, clust, titletxt = "", pal, crs = albers_crs) {

  stopifnot(nrow(coords) == length(clust))
  
  # Build sf points; name the cluster column and make it a factor
  pts <- data.frame(x = coords$X, y = coords$Y, cluster = factor(clust))
  if (!is.null(names(pal)) && length(names(pal))) {
    pts$cluster <- factor(pts$cluster, levels = names(pal))  # match palette names
  }
  sf_pts <- st_as_sf(pts, coords = c("x", "y"), crs = crs)
  
  ggplot() +
    geom_sf(data = sf_pts, aes(color = cluster), size = 0.9) +
    geom_sf(data = coast, fill = NA, color = "grey10", linewidth = 0.5, inherit.aes = FALSE) +
    scale_color_manual(values = pal, limits = levels(sf_pts$cluster), drop = FALSE, name = "Cluster") +
    labs(title = titletxt, x = NULL, y = NULL) +
    coord_sf(expand = FALSE) +
    theme_minimal() +
    theme( # theme grows boxes and text but not glyph size
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))
}


#---- MakeScreePlot: returns a ggplot. ----
# samp is optional, uses all dat if omitted.
MakeScreePlot <- function( indat, nclust, nrand, maxi, sampsize = 0 ){
  #initialize list for results
  wss <- numeric(nclust) 
  
  #subsample as requested
  if (sampsize > 0) {
    samp <- sample( 1:length( indat[ , 1] ), sampsize )
    dat <- indat[ samp, ]
  } else dat <- indat
  
  for (i in 1:nclust) {
    # Fit the model: km.out
    print( paste0( "centers ",i))
    km.out <- kmeans(dat, centers = i, nstart = nrand, iter.max = maxi)
    # Save the within cluster sum of squares
    wss[i] <- km.out$tot.withinss
    
    # calculate the silhouete width
  }
  
  # Produce the scree plot ... using a tibble, I guess. 
  wss_df <- tibble(clusters = 1:nclust, wss = wss)
  scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    xlab('Number of clusters') +
    ylab('Total within-cluster sum of squares') +
    geom_hline(
      yintercept = wss, 
      linetype = 'dashed')
  
  return( scree_plot)
}


#---- ClusterPCA: Returns a pair of PCA plots showing the separation of the clusters.
# NOTE: Uses a relatively small subset of the overall data so these will change with a different sample.
ClusterPCAall <- function( p_data ) {
  
 # ssidx <- sample( 1:length( stack_data_clean[ , 1] ), n_samp )
 # ssamp <- stack_data_clean[ ssidx, ]
  
  # re-run cluster for smaller sample.
  #cluster_result <- kmeans(datssamp, centers = clustnum, nstart = randomz) # less than 10 seconds
  #csamp <- cluster_result$cluster
  
  # Create the data structure for PCA profiling (combining predictor data with the clusters). 
  # Put cluster # first so easy to find for PCA
  #p_data <- as.data.frame( cbind(cluster = clust, dat ) )
  
  res_pca <- prcomp(p_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
  ind_coord <- as.data.frame(get_pca_ind(res_pca)$coord)
  # Add clusters from the classification
  ind_coord$cluster <- factor( p_data$cluster )
  # Data inspection
  #head(ind_coord)
  
  # Percentage of variance explained by dimensions
  eigenvalue <- round(get_eigenvalue(res_pca), 1)
  var_percent <- eigenvalue$variance.percent
  
  # Look at the clusters for the first 4 PCs
  a <- ggscatter(
    ind_coord, x = "Dim.1", y = "Dim.2", 
    color = "cluster", palette = our_pal, ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 1 (", var_percent[1], "% )" ),
    ylab = paste0("Dim 2 (", var_percent[2], "% )" )
  ) +
    stat_mean( aes(fill=cluster), shape=21, color="black", size = 4)

  pca_loadings <- data.frame(res_pca$rotation)
  pca_scores   <- data.frame(res_pca$x)
  pca_loadingsA <- pca_loadings * max(abs(pca_scores$PC1), abs(pca_scores$PC2))

  plot1 <- a +
    geom_segment(data = pca_loadingsA, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = pca_loadingsA, aes(x = PC1, y = PC2, label = rownames(pca_loadingsA)), 
              hjust = 0, vjust = 1, color = "black")

  a <- ggscatter(
    ind_coord, x = "Dim.3", y = "Dim.4", 
    color = "cluster", palette = our_pal, ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 3 (", var_percent[3], "% )" ),
    ylab = paste0("Dim 4 (", var_percent[4], "% )" )
  ) +
    stat_mean( aes(fill=cluster), shape=21, color="black", size = 4)
  
  pca_loadingsB <- pca_loadings * max(abs(pca_scores$PC3), abs(pca_scores$PC4))
  
  plot2 <- a +
    geom_segment(data = pca_loadingsB, aes(x = 0, y = 0, xend = PC3, yend = PC4),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = pca_loadingsB, aes(x = PC3, y = PC4, label = rownames(pca_loadingsB)), 
              hjust = 0, vjust = 1, color = "black")
  
  print( "Done clusters ...")
  return( list(res_pca, plot1, plot2) )
}

#---- ClusterPCA: Returns a pair of PCA plots showing the separation of the clusters.
# NOTE: Uses a relatively small subset of the overall data so these will change with a different sample.
ClusterPCA <- function( n_samp, clustnum ) {
  
  ssidx <- sample( 1:length( stack_data_clean[ , 1] ), n_samp )
  ssamp <- stack_data_clean[ ssidx, ]
  
  # re-run cluster for smaller sample.
  cluster_result <- kmeans(ssamp, centers = clustnum, nstart = randomz) # less than 10 seconds
  csamp <- cluster_result$cluster
  
  # Create the data structure for PCA profiling (combining predictor data with the clusters). 
  # Put cluster # first so easy to find for PCA
  p_data <- as.data.frame( cbind(cluster = csamp, ssamp ) )
  
  res_pca <- prcomp(p_data[,-1],  scale = TRUE)
  # PC coordinates of individual raster cells
  ind_coord <- as.data.frame(get_pca_ind(res_pca)$coord)
  # Add clusters from the classification
  ind_coord$cluster <- factor( p_data$cluster )
  # Data inspection
  #head(ind_coord)
  
  # Percentage of variance explained by dimensions
  eigenvalue <- round(get_eigenvalue(res_pca), 1)
  var_percent <- eigenvalue$variance.percent
  
  # Look at the clusters for the first 4 PCs
  a <- ggscatter(
    ind_coord, x = "Dim.1", y = "Dim.2", 
    color = "cluster", palette = our_pal, ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 1 (", var_percent[1], "% )" ),
    ylab = paste0("Dim 2 (", var_percent[2], "% )" )
  ) +
    stat_mean( aes(fill=cluster), shape=21, color="black", size = 4)
  
  pca_loadings <- data.frame(res_pca$rotation)
  pca_scores   <- data.frame(res_pca$x)
  pca_loadingsA <- pca_loadings * max(abs(pca_scores$PC1), abs(pca_scores$PC2))
  
  plot1 <- a +
    geom_segment(data = pca_loadingsA, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = pca_loadingsA, aes(x = PC1, y = PC2, label = rownames(pca_loadingsA)), 
              hjust = 0, vjust = 1, color = "black")
  
  a <- ggscatter(
    ind_coord, x = "Dim.3", y = "Dim.4", 
    color = "cluster", palette = our_pal, ellipse = TRUE, ellipse.type = "convex",
    size = 1.5,  legend = "right", ggtheme = theme_bw(),
    xlab = paste0("Dim 3 (", var_percent[3], "% )" ),
    ylab = paste0("Dim 4 (", var_percent[4], "% )" )
  ) +
    stat_mean( aes(fill=cluster), shape=21, color="black", size = 4)
  
  pca_loadingsB <- pca_loadings * max(abs(pca_scores$PC3), abs(pca_scores$PC4))
  
  plot2 <- a +
    geom_segment(data = pca_loadingsB, aes(x = 0, y = 0, xend = PC3, yend = PC4),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = pca_loadingsB, aes(x = PC3, y = PC4, label = rownames(pca_loadingsB)), 
              hjust = 0, vjust = 1, color = "black")
  
  print( "Done clusters ...")
  return( list(res_pca, plot1, plot2) )
}

#---- PredictClusters: returns cluster assignments for un-clustered pictures. ----
# Assembly required with the classified pixels before plotting. Function by ChatGPT.
PredictClusters <- function(newdata, kmeans_model) {
  centers <- kmeans_model$centers
  dists <- as.matrix(dist(rbind(centers, newdata)))
  dists <- dists[-(1:nrow(centers)), 1:nrow(centers)]
  cluster_assignments <- apply(dists, 1, which.min)
  return(cluster_assignments)
}


#---- TransferClusters: returns new values for the prediction rasters containing predictions. ----
# transferred from the random sample used in the classification ----
# Assembly required with the classified pixels before plotting. Function by ChatGPT.
transferCluster <- function(values_target, c_result){
# uses GLOBALS sidx, stack_data_clean   
  #--- Two steps here: First assign clusters to the rest of the clean stack data, 
  #     THEN put the clean data back in the raster.
  # Uses samp and sidx from lines ~150 above.
  # Find the things not in the sample (using sidx from line ~150 above).
  
  not_sidx <- setdiff( 1:dim(stack_data_clean)[[1]], sidx )
  
  # the data to predict clusters for
  new_dat <- stack_data_clean[ not_sidx, ]
  # Process the new data in chunks of 10k to ensure performance
  chunk_size <- 10000
  n_chunks <- ceiling(nrow(new_dat) / chunk_size)
  
  print( paste0( "predicting clusters for ", dim(new_dat)[[1]], " pixels using ", 
                 n_chunks, " chunks of ", chunk_size, ". Stand by ... ") )
  
  # Where we are putting our new chunks
  predicted_clusters <- vector("integer", nrow(new_dat))
  
  for ( i in seq_len(n_chunks) ) {
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, nrow(new_dat))
    
    chunk <- new_dat[start_index:end_index, ]
    predicted_clusters[start_index:end_index] <- PredictClusters(chunk, c_result)
    cat(i, " ")
  }
  
  #length(predicted_clusters)
  #length(c_result$cluster)
  # Combine the results for the cleaned data ... 
  clean_clusts <- array( 1:dim(stack_data_clean)[1] )
  clean_clusts[ sidx ] <- c_result$cluster
  clean_clusts[ not_sidx ] <- predicted_clusters
  
  # replace clustered values with the cluster results
  values_target[ clean_idx ] <- clean_clusts
  # reassign to the plotting raster
  return( values_target )
}





####------------Depreciated or replaced functions----------------------------

#---- TrimStack: limit extents based on a polygon mask (vestigial)
# NOTE: This function now superseded by the crop/mask approach in ClipPredictors.
TrimStack <- function( stack_in, padsize ) {
  # Assigning extents here worked but results odd - plots appeared but looked empty.
  y <- stack()
  # estimate extents from first raster ... 
  x <- trim(stack_in[[1]], padsize )
  x_ext <- extent(x)
  x_ext <- ceiling( x_ext )
  extent( x ) <- x_ext
  y <- stack(y, x)
  print( paste0('layer 1 trimmed.'))
  
  # use these extents for the remaining layers ...
  for (i in 2:dim( stack_in)[[3]] ) {
    x <- setExtent( stack_in[[i]], x_ext, keepres=T, snap=T )
    y <- stack( y, x )
    print( paste0('layer ', i, ' trimmed.'))
  }
  return( y )  
}


# if (spacesub) {
#   
#   sMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
#   smMask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
#   
#   # Data extents No mask and polygons
#   x  <- prepped_layers
#   
#   # Use only one of the loaded shapefiles
#   x  <- raster::mask(x, sMask )
#   # x  <- raster::mask(x, smMask )
#   
#   # Plot first layer in the resulting stack for inspection
#   plot(trim(x[[1]]), main = 'Local')
#   
#   # Plot masks over raster. NB: will only see if raster extends beyond the masks. 
#   sp::plot( sMask, border = "darkgreen", lwd = 2, add=TRUE)
#   sp::plot( smMask, border = "red", lwd = 2, add=TRUE)
#   
#   prepped_layers <- x
#   rm('x')
# }



ScalePredictors <- function( the_stack ){
  # A little shimmy to avoid scaling substrate and name it consistently
  #  scaled_stack <- scale( dropLayer( the_stack, "SUBSTRATE") )
  #  scaled_stack <- stack( scaled_stack, trim_layers$SUBSTRATE )
  # safe to assume its the last layer in the stack
  #  names( scaled_stack[[ dim(scaled_stack)[[3]] ]] ) <- "substrate"
  
} 


# Apply scaling to standardize rasters. Does either min/max or z-score.
# Returns a scaled raster.
Scale.Raster <- function( in.raster, scale.how = "mm", intshift = 'T' ) {
  
  if (scale.how == "z") {
    mean.value <- cellStats( in.raster, stat='mean', na.rm=TRUE )
    sd.value   <- cellStats( in.raster, stat='sd', na.rm=TRUE )  
  } else {
    old.min    <- cellStats( in.raster, stat='min', na.rm=TRUE )  
    old.max    <- cellStats( in.raster, stat='max', na.rm=TRUE )  
  }
  
  if (scale.how == "z") {
    # Perform z-score normalization
    scaled.raster <- (in.raster - mean.value) / sd.value
  } else {
    # Perform min/max scaling
    scaled.raster <- (in.raster - old.min) / (old.max - old.min)
  }
  
  if (intshift) {
    scaled.raster <- as.integer( scaled.raster * 1000 )
  }
  
  return( scaled.raster )
}


# fin
