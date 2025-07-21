#################################################################################
# Script:  Classification_main_LSSM.R - Broughton LSSM version
# Created: February 2024. EJG
# 
# This script sources the necessary libraries and functions, coordinates the 
# analysis, and creates data structures that are then 'knitted' together (I think).
# So the idea is to run bits in here, and then 'Render' the RMD script. 
# Seems straightforward. :)
#
## Updates: 
# 2024/04/29: Steady development past few weeks; alll necessary pieces now developed. 
# 2024/04/29: Git repo created and working version pushed prior to RMarkdown development.
# 2024/05/02: Completed smoothing pass thru code; sorted raster plotting. Pushed.
# 2024/05/07: Another pass thru, adding some controls. Ready for RMD work.Pushed.
# 2024/05/29: After some exploratory work with Romina's tifs from SPECTRAL lab, forked the 
#   code so this version focuses on the LSSM needs.
# 2024/06/05: Updated with relevant code from latest DFO version.
# 2024/08/12: Pre-holidays, worked thru code and began work on rMarkdown script. This script
#   is now ahead of the DFO version. Refamiliarization in progress. 
# 2024/08/30: Back to work on classification. Began by renaming the 2 projects for clarity, and 
#   moved all progress over to DFO version -> finish that first.
# 2024/09/22: DFO work close to done. Nicely developing RMD. Back to here to create 
#  initial clusters for the LSSM. 
#  - 2 key bits to integrate: distributions and scaling, and RMD report.
# 2024/09/22: Lots of RMD work. Acceptable PDF report of cluster methods/results.
#  - GIT Repos and scripts renamed for clarity.
# 2024/10/02: Further updates to RMD. Can now run reasonable cluster reports.
#  - exercising thoroughness by examining full basket of predictors ... 
# 2024/10/04: Formal walk-thru to winnow predictors.
# 2025/07/15: Start work to complete with original Biannucci data and learnings from DFO verson.


# TO DO: 
#  1. Reproduce (ish) Romina - use as prototype
#  2. Combine attributes and develop update plan
#  4. Spend a TINY bit of time to see if MSEA data can "improve" fit.
#  5. Can repeated clusterings of same data be attached to same RMD report,
#     OR easier by hand, post-hoc?
#######################################################################################

print('Starting Classification  ... LSSM Version')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "classification_functions.R" )

today <- format(Sys.Date(), "%Y-%m-%d")

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
raster_dir <- 'C:/Data/SpaceData/Classification/Romina'
data_dir   <- 'C:/Data/Git/Classification-LSSM/Data'
results_dir<- 'C:/Data/Git/Classification-LSSM/Results' 

# Processing FLAGS...
loadtifs <- F # If true the data will be re-loaded from TIFs, else it will be loaded from rData.
clipdata <- T # If true a spatial subset of the data will be taken based on a polygon shape file. 
reclust  <- T # If true, re-cluster full data set prior to mapping, else predict to unclassified pixels.

#---- Part 1: Load and trim predictor data.  ----
# If loadtifs == TRUE then TIFs loaded from raster_dir, else a stack is loaded from rData file.

tif_stack <- stack()

if (loadtifs) {
  print( "Loading predictors ... ")
  src_stack <- LoadPredictors( raster_dir )
  print( "Data loaded.")
  
  tif_stack <- src_stack
  
  if (clipdata) {
    print( "clipping TIFs ... ")
    #amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_small.shp")
    amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smallest.shp")
    tif_stack <- ClipSPECPredictors( tif_stack, amask )
    print('Rasters clipped.')
  }
 
  
  writeRaster( src_stack[[3]], paste0( results_dir, "/tif_test.tif"), overwrite=TRUE)
  
# IF we want to include DFO data it will need to match Bianucci's spatial reference.

  save( tif_stack, file = paste0( data_dir, '/tifs_SPECTRAL_clipped_', today, '.rData' ))
  
} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
  load( paste0( data_dir, '/tifs_SPECTRAL_clipped_2025-07-14.rData' ))
}


#---- Part 2: Predictor correlations ----
#-- Move to matrix space from raster space 
x <- getValues( tif_stack )

#-- Identify and use only pixels with all data
# NB: This decouples the data from the RasterStack, and requires re-assembly for mapping
dim(x)
clean_idx <- complete.cases( x )
x_clean <- x[ clean_idx, ]
dim( x_clean )

#---- Correlation across data layers 
cor_table <- round( cor( x_clean ), 3)
cor_table[lower.tri(cor_table, diag=TRUE)] <- NA

high_rows <- apply(cor_table, 1, function(row) any(row > 0.6, na.rm = TRUE))
z <- cor_table[ high_rows, ]


tx_clean <- x_clean


stack_data <- getValues( tif_stack)
print( "Transform the data  ... ")
# Transform those predictors with abs(skew) > 1. Predictor selection was done manually and 
# includes 10 of 17: 
#   julBT_min, julBS_min, julSS_min, julSS_max,
#   julBSpd_min, julBSpd_max, julSSpd_min, julSSpd_max,
#   roughness, tidal_cur

# MakeMoreNormal() is hard-coded based on trial and error. 
# trans parameter not used but could be in order to easily include the transform applied 
# to the resulting table, but do the work only if necessary. 
t_stack <- MakeMoreNormal(stack_data, trans )
colnames(t_stack)


# Scale and center the transformed data.
print("Centering and scaling  ... ")
tmp_stack <- scale(t_stack, center = T,  scale = T)
t_stack_data <- as.data.frame(tmp_stack)
print('Data centered.')


# Sequential dropping of predictors based on analytic results. 
# Done here so that correlation plots are updated.  
# Step 1: Remaining correlated predictors from the SPECTRAL data set dropped.
# Step 2: Drop Northness as it continues to be uniformly distributed across clusters.
# Step 3: Drop tidal_cur as performance continues to be highly correlated with julSSpd_min
t_stack_data <- 
  t_stack_data[, !colnames(t_stack_data) %in% 
      c("julBT_min","julST_min","julSSpd_min","julBSpd_max", "julBS_max", "julBS_min", "julSS_min",
        "northness", "tidal_cur") ]

colnames( t_stack_data )

# A quick correlation across data layers
x <- t_stack_data
x_clean <- x[ complete.cases(x), ]
cor_table <- cor( x_clean )
cor_table[lower.tri(cor_table)] <- NA
b_cor <- (abs(cor_table) >= 0.6) & (cor_table != 1)
sum( b_cor, na.rm=T)

# remove any rows with an NA
# NB: This decouples the data from the RasterStack and requires re-assembly for mapping
# THESE are the two key data structures used in subsequent steps

clean_idx <- complete.cases(t_stack_data)
stack_data_clean <- t_stack_data[ clean_idx, ]

dim( t_stack_data )
dim( stack_data_clean )



#---- Part 2 of 3: Cluster number selection ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

#---- Part 2a: Explore number of clusters using Within-sum-of-squares scree plot ----
# Runs kmeans with increasing number of clusters

nclust   <- 18 # number of clusters for scree plot
nsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( stack_data_clean, nclust, randomz, imax, nsample )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 6 # the number of clusters based on scree plot, above.
nsample <- 500000 # a larger sample for more robust classification

sidx <- sample( 1:length( stack_data_clean[ , 1] ), nsample )
samp <- stack_data_clean[ sidx, ]


## JUNE 2025 - working here ... 
# Run all the data ... Can do for KKD.

samp <- stack_data_clean

cluster_result <- kmeans(samp, centers=nclust, nstart=randomz, iter.max=imax) 

#---- Part 2b: Create heat map of within-cluster standard deviations ----

# Define color palette
pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette

profile_data <- as.data.frame( cbind(cluster = cluster_result$cluster, samp ) )

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
head(x)
xm <- melt( x, id.var = "cluster" )

z_heat <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "", x = "Clusters", y = "Predictors", fill = "Value")
z_heat

#---- Part 2c: Examine silhouette plot of the WORKING clusters  ----
# Uses the predictor values and the corresponding assigned cluster
# Need to subsample from the cluster result above as distance matrix take long time.

# Take a subsample of the clusters and the predictors for the silhouette plot. 
sil_n <- 20000
silx <- sample( 1:length( samp[ , 1] ), sil_n )

cs <- cluster_result$cluster[ silx ]
ss <- samp[ silx, ]

# Calculate a distance matrix between the predictors and plot assigned to each cluster.
# both these steps are time consuming, hence a smaller sample.
c_dist <- dist(ss)
sk <- silhouette(cs, c_dist)
#mean( sk[,"sil_width"] )

plot(sk, col = 1:nclust, border=NA, main = "" )


#---- Part 3: Detailed examination of N clusters  ----
#---- Part 3a: Show cluster groupings using PCA ----

#-- Can take some time so it makes its own cluster 
pca_n <- 25000

pca_results <- ClusterPCA( pca_n, nclust ) # uses global variable stack_data_clean
names(pca_results) <- c("loadings", "plot1","plot2")

pca_results$plot1
#Percentage of variance explained by dimensions
#eigenvalue <- round(get_eigenvalue(res_pca), 1)
#var_percent <- eigenvalue$variance.percent

#---- Part 3b: Violins of predictor contributions to WORKING clusters ----

x <- as.data.frame( samp )
x$cluster <- as.factor( cluster_result$cluster )

y <- x %>%
  pivot_longer(cols = -cluster, names_to = "predictor", values_to = "value")

# Create violin plot
vplots <- 
  ggplot(y, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ predictor, scales = "free_y") +
  theme_minimal() +
  labs(title = "Violin Plots of Predictors Across k-means Clusters",
       x = "Cluster",
       y = "Value")
vplots

#---- Part 3c: Spatialize and display the WORKING clusters ----
# NB: To show a comprehensive map, can either:
#     a) re-cluster the entire data set (using imax and randomz from above) or
#     b) Predict to the unsampled portion of the raster. 


# For the smallest extents, don't need to recluster (below)

cluster_raster <- tif_stack[[1]]
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

#   # Assign the clustered values ... 
#   # extract values from the target cluster
   new_values <- values( cluster_raster )
#   # replace non-NA values with the cluster results
   new_values[ clean_idx ] <- cluster_result$cluster
#   # put the updated values back on the target cluster
   values( cluster_raster ) <- new_values  

writeRaster( cluster_raster, paste0( results_dir, "/tif_stack_shallow_", today, ".tif"), overwrite=TRUE)


### Jun 25 '25: Reconstruction didn't work with the 'smallest' clip of the data, 
# likely cuz it can all be classified so this logic doesn't quite work. 
# NEED TO build a condition in so that its only used when needed. 

## initialize target data structure 
# cluster_raster <- tif_stack[[1]]
# dataType(cluster_raster) <- "INT1U"
# cluster_raster[] <- NA
# 
# if (reclust == T) {
#   # Re-cluster using all the clean data  ... 
#   # less than 1 min with iter.max = 20, nstart = 20 for smallest region
#   
#   # Re-use seed from above to ensure identical clusters
#   .Random.seed <- saved_seed
#   cluster_result <- kmeans(stack_data_clean, centers = nclust, nstart = randomz, iter.max = imax)
#   
#   # Assign the clustered values ... 
#   # extract values from the target cluster
#   new_values <- values( cluster_raster )
#   # replace non-NA values with the cluster results
#   new_values[ clean_idx ] <- cluster_result$cluster
#   # put the updated values back on the target cluster
#   values( cluster_raster ) <- new_values  
# } else {
#   # Predict values for unclustered cells. 
#   # Currently more time-consuming than re-classifying everything, probably cuz predicting
#   # includes all the pixels outside the area of interest. Could likely be more efficient. 
#   values( cluster_raster ) <- transferCluster( tx_clean, values(cluster_raster), cluster_result )
# }


#--- Display the results, first as histogram then as map.  
raster::hist( values(cluster_raster ))

# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent") # Max for Accent is 8

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
z_map <- levelplot( cluster_raster, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )
z_map
writeRaster( cluster_raster, paste0( results_dir, "/SPECtrim_testing_", today, ".tif"), overwrite=TRUE)


#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
rmarkdown::render( "Classification_LSSM_PDF.Rmd",   
#rmarkdown::render( "Classification_test.Rmd",   
                   output_format = 'pdf_document',
                   output_dir    = results_dir,
                   output_file = paste0( "LSSM_SPECtrim_testing_", today ))


#---- Some details on correlation analysis ... ----
#-- Correlation across UN-scaled data layers ... 
# foo <- getValues( scaled_layers )
# foo_clean <- na.omit(stack_data)
# pick <- sample( 1:length( foo_clean[ , 1] ), 10000 )
# plot( foo_clean[ pick,'rugosity'] ~ foo_clean[ pick,'standard_deviation_slope'] )

#---- ----

# FIN.






