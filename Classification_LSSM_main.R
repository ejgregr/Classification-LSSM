# Script:  Classification_main_LSSM.R - Broughton LSSM version
# Created: February 2024. EJG
############################################################################----
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
# 2025/07/15: Start work to complete with original Bianucci data and learnings from DFO version.
# 2025/07/21: Completed processing of original Bianucci data to points
#   --> Start major revision to classification, replacing Spectral rasters with DFO FVCOM points.
# 2025/07/22: Completed first pass thru LSSM classes. Pushed. 
# 2025/09/03: Back at this now, tidying up loose ends, and coming up with a plan.


# To DO:
# - Examine predictors with cor = 1.0. FVCOM or data transfer problem?

################################################################################

print('Starting Classification  ... LSSM Version 2: Points not rasters')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
 source( "classification_functions.R" )
# Only if you need to reload FVCOM data 
#source( "netCDF_processing.R" )

today <- format(Sys.Date(), "%Y-%m-%d")

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
source_dir <- 'C:/Data/SpaceData/Broughton/netcdf'
data_dir   <- 'C:/Data/Git/Classification-LSSM/Data'
results_dir<- 'C:/Data/Git/Classification-LSSM/Results' 

# NOTES: 
# - For loading of TIFs, and the necessary subsampling, see the DFO version of the classification code.
# - The FVCOM data are only sub-sampled spatially, in ArcGIS.

#---- Load FVCOM point data ----
# -Jul 21 version of FVCOM clipped to KD extents
#load( file = paste0( source_dir, '/FVCOM_point_data_2025-07-22.rData' ))

# Sept 5 version. Minor code tweaks
load( file = paste0( data_dir, '/FVCOM_point_data_2025-09-05.rData' ))

str(kd_pts)
fv_dat <- kd_pts

# Shorten some names ... 
names(fv_dat) <- gsub("salinity", "salt", names(fv_dat))
names(fv_dat) <- gsub("bottom",   "bott", names(fv_dat))
names(fv_dat) <- gsub("surface",  "surf", names(fv_dat))
colnames(fv_dat)

# remove (X,Y) for the analysis
#plot( kd_pts$X, kd_pts$Y)
fv_dat <- fv_dat[ , !(names(fv_dat) %in% c("X", "Y")) ]

#---- Predictor correlations ----
cor_table <- round( cor( fv_dat ), 3)
cor_table[lower.tri(cor_table, diag=TRUE)] <- NA

high_rows <- apply(cor_table, 1, function(row) any(row > 0.5, na.rm = TRUE))
z <- cor_table[ high_rows, ]
showHighestPairs2( z, 0.5 )

# Start removing variables. 

# Pass 1:
# These are perfectly correlated. Odd?
# "salt_DEC_bott_ave" "temp_DEC_bott_ave" "salt_DEC_surf_max"
# "salt_DEC_surf_ave" "temp_DEC_surf_ave" "CS_DEC_surf_ave"  
# plot( fv_dat$temp_DEC_surf_ave, fv_dat$temp_DEC_bott_ave)
# Drop bottoms.
fv_dat <- fv_dat[ , !grepl("salt_DEC_bott_ave", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("temp_DEC_bott_ave", names(fv_dat)) ]
# Retain chemistry over currents.
fv_dat <- fv_dat[ , !grepl("CS_DEC_surf_ave", names(fv_dat)) ]

# Pass 2: 
# 15 correlations >= 0.9. First drop correlated averages
fv_dat <- fv_dat[ , !grepl("temp_DEC_surf_ave", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("salt_DEC_surf_ave", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("salt_JUL_bott_ave", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("temp_JUL_bott_ave", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("temp_JUL_surf_ave", names(fv_dat)) ]
# min correlated with max, prefer max
fv_dat <- fv_dat[ , !grepl("temp_JUL_surf_min", names(fv_dat)) ] 
# remainder all correlated with Dec CS. Retain surface where possible 
fv_dat <- fv_dat[ , !grepl("CS_DEC_bott_min", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("CS_DEC_bott_max", names(fv_dat)) ]
# Dec/Jan surf current max and mins correlated. Keep JUL
fv_dat <- fv_dat[ , !grepl("CS_DEC_surf_min", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("CS_DEC_surf_max", names(fv_dat)) ]
# Retain min over average
fv_dat <- fv_dat[ , !grepl("salt_JUL_surf_ave", names(fv_dat)) ]

#Pass 3:
# 8 correlations >= 0.8. First drop correlated averages
#  retain max on temp, and min on salt ...   
fv_dat <- fv_dat[ , !grepl("temp_DEC_bott_min", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("salt_JUL_bott_max", names(fv_dat)) ]
# Remaining 6 are all JUL currents 
fv_dat <- fv_dat[ , !grepl("CS_JUL_surf_min", names(fv_dat)) ] # correlated w 3
fv_dat <- fv_dat[ , !grepl("CS_JUL_bott_min", names(fv_dat)) ]
# check if either further cor'd but not so drop bottom
fv_dat <- fv_dat[ , !grepl("CS_JUL_bott_max", names(fv_dat)) ]

#Pass4:
# 9 correlations >= 0.6. Most min/max pairs, drop min or max depending
fv_dat <- fv_dat[ , !grepl("salt_DEC_bott_max", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("temp_DEC_surf_min", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("temp_JUL_bott_min", names(fv_dat)) ]
# keep a bottom value for DEC
fv_dat <- fv_dat[ , !grepl("temp_DEC_surf_max", names(fv_dat)) ]
# min salts at bottom and surface correlated across seasons. 
# retain Jul surface and Dec bot
fv_dat <- fv_dat[ , !grepl("salt_JUL_bott_min", names(fv_dat)) ]
fv_dat <- fv_dat[ , !grepl("salt_DEC_surf_min", names(fv_dat)) ]

names(fv_dat)
dim(fv_dat)

hist1 <- PlotHistos( fv_dat, "Untransformed, uncorrelated FVCOM predictors" )
hist1
# Above correlation analysis done using untransformed predictors. 
# Now work on Normality 

#---- Transformation and scaling of uncorrelated predictors  ----
print( "Transform the data  ... ")

# Show predictors with abs(skew) > 1.
skews <- apply(fv_dat, 2, skewness, na.rm = TRUE)
skews[skews > 1]

# MakeMoreNormal() is hard-coded based on trial and error. 
# Adapted to FVCOM data
tfv_dat <- MakeMoreNormal( fv_dat )

dim(tfv_dat)
hist2 <- PlotHistos( tfv_dat, "Transformed FVCOM predictors" )
hist2

# Now scale and center the transformed data.
print("Centering and scaling  ... ")
x <- scale( tfv_dat, center = T,  scale = T )
done_dat <- as.data.frame( x )
hist3 <- PlotHistos( done_dat, "Transformed and scaled FVCOM predictors" )
hist3

#---- Final data selection ---- 
# Any sequential dropping of predictors based on analytic results done here
# So all clustering reflects the changes. 

use_dat <- done_dat

# Start with a seasonal comparison

# Drop JuL for winter data only - yields only 3 predictors.
use_dat <- use_dat[ , !grepl("JUL", names(use_dat)) ]

# Drop current data to focus on water chemistry
# use_dat <- use_dat[ , !grepl("CS", names(use_dat)) ]

# * This is the df used throughout below *
names( use_dat )


#---- Cluster analysis Part 1 - Number of clusters  ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

# Explore number of clusters using Within-sum-of-squares scree plot
# Runs kmeans with increasing number of clusters
nclust     <- 18 # number of clusters for scree plotnsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
scree_plot <- MakeScreePlot( use_dat, nclust, randomz, imax, 0 )
scree_plot


#---- Cluster analysis Part 2 - Create working set of N clusters ----

#-- SETUP for clustering and figure labeling ----
# confirm data subset as desired 
use_dat <- done_dat[ , !grepl("DEC", names(done_dat)) ]
nclust  <- 9 # the number of clusters based on scree plot, above.

clust_note <- "seasonal_JUL_"
title_text <- "9 clusters"

cluster_result <- kmeans(use_dat, centers=nclust, nstart=randomz, iter.max=imax) 

# Setup text for plot and file names 

# Save the current cluster.
JUL_9clust <- cluster_result


#---- Part 3: Cluster analysis diagnostic plots  ----

# The cluster and data to use in the plots
use_clust <- JUL_9clust
use_dat   <- use_dat

#---- Heat map of within-cluster standard deviations ----
heat_map <- PlotHeatMap( use_clust$cluster, use_dat, title_text  ) 
heat_map
#---- Examine silhouette plot of the WORKING clusters  ----
c_dist <- dist( use_dat )
sk <- silhouette(use_clust$cluster, c_dist)
sil_plot <- plot(sk, col = 1:nclust, border=NA, main = "" )

#---- Show cluster groupings using PCA ----
pca_results <- ClusterPCAall( cbind(cluster = use_clust$cluster, use_dat ) )
names( pca_results ) <- c("loadings", "plot1","plot2")

#Percentage of variance explained by dimensions
eigenvalues <- round(get_eigenvalue(pca_results[[1]]), 2)
var_percent <- eigenvalues$variance.percent

#---- Violins of predictor contributions to WORKING clusters ----
v_plots <- PlotViolins( use_clust$cluster, use_dat )


#---- Part 4: Map the working point clusters ----
# Jul22: Currently using points here. 
# For a raster, see the code in the DFO classification version.

pts <- data.frame(x = kd_pts$X, y = kd_pts$Y, cluster = use_clust$cluster)
# 32609 is the CRS for UTM Zone 9N, the projection of the FVCOM data
sf_pts <- st_as_sf(pts, coords = c("x", "y"), crs = 32609)

theme_set(theme_gray()) # Reset ggplot() defaults
clust_map <- ggplot(sf_pts, aes(color = factor(cluster))) +
  geom_sf( size=0.8) + theme_minimal() +
  labs(color = "Clusters") +
  guides(color = guide_legend(override.aes = list(size = 2))) +  # controls legend point size
  theme(legend.position = "bottom") +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.key.size = unit(0.6, "lines"),    # smaller legend keys (boxes/symbols)
    legend.text = element_text(size = 8),    # smaller legend labels/text
    legend.title = element_text(size = 10) )   # (optional) smaller legend title)
#  guides(color = guide_legend(override.aes = list(size = 2))) +  # controls legend point size
#  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
#  
clust_map

# Output a shape file
# sf_pts created above
outname <- paste0( "LSSM_Results_v", as.character(dim(use_dat)[[2]]), "_c", nclust )
st_write(sf_pts, paste0( results_dir, outname, ".shp"), delete_layer = TRUE)

# Define color palette - Max N for Accent is 8
#pal_clust <- brewer.pal(n = nclust, "Accent") 




#---- Standardize pairs of clusters ----


km_ref <- JUL_9_clust
km_new <- DEC_9_clust

# Cross-tab BEFORE alignment
cat("\nContingency (ref vs new) BEFORE alignment:\n")
print( table(km_ref$cluster, km_new$cluster) )

# ---- Align by membership overlap ----
align <- match_clusters_by_overlap(km_ref$cluster, km_new$cluster)

cat("\nMapping (index = new label, value = ref label):\n")
print(align$mapping)

cat("\nContingency used for assignment:\n")
print(align$contingency)

# Check AFTER alignment
cat("\nContingency (ref vs aligned) AFTER alignment:\n")
print(table(km_ref$cluster, align$aligned))



# Prep the FVCOM data ... 
names( done_dat )
str(sf_pts)

km_ref <- JUL_9_clust
km_new <- DEC_9_clust

# Standardize the predictor names ... 
colnames(km_ref$centers) <- gsub("_JUL_", "_", colnames(km_ref$centers), fixed = TRUE)
colnames(km_new$centers) <- gsub("_DEC_", "_", colnames(km_new$centers), fixed = TRUE)

# ---- Usage ----
# Assume you already have:
#   km_ref   <- kmeans(X_ref, centers = K, nstart = 25)
#   km_new   <- kmeans(X_new, centers = K, nstart = 25)
# and a coords data frame `pts` with columns x,y in the SAME row order as X_ref/X_new.

aln <- align_kmeans_by_overlap(km_ref, km_new)

# Check alignment quality
table(km_ref$cluster, aln$aligned)


# Random, vivid, distinct colors (HCL space)
random_pal <- function(K, seed = NULL, l = 65, c = 100) {
  if (!is.null(seed)) set.seed(seed)
  hues <- sample(seq(0, 360 - 360 / K, length.out = K))
  stats::setNames(grDevices::hcl(h = hues, c = c, l = l), as.character(1:K))
}

# Example:
pal <- random_pal(K, seed = 42)   # set seed for reproducible "random"

# pts must be the same row order used for kmeans
p_ref <- ggplot(pts, aes(x, y, color = factor(km_ref$cluster))) +
  geom_point(size = 1) + scale_color_manual(values = pal, drop = FALSE) +
  coord_equal() + ggtitle("Reference")

p_new <- ggplot(pts, aes(x, y, color = factor(aln$aligned))) +
  geom_point(size = 1) + scale_color_manual(values = pal, drop = FALSE) +
  coord_equal() + ggtitle("New (aligned to reference)")

print(p_ref); print(p_new)

# Or side-by-side if you have {patchwork}:
# (p_ref | p_new)








#---- Output images as png ----
outname <- paste0( results_dir, "/", clust_note, "_scree_plot_v", as.character(dim(use_dat)[[2]]), ".png" )
png(outname, width = 800, height = 600, res = 200)  # open PNG device
scree_plot
dev.off()    

outname <- paste0( results_dir, "/", clust_note, "_heat_map_v", as.character(dim(use_dat)[[2]]), "_c", nclust, ".png" )
png(outname, width = 800, height = 600, res = 200)  # open PNG device
heat_map
dev.off()    

outname <- paste0( results_dir, "/", clust_note, "_violin_plots_v", as.character(dim(use_dat)[[2]]), "_c", nclust, ".png" )
png(outname, width = 1200, height = 800, res = 200)  # open PNG device
v_plots
dev.off()    

outname <- paste0( results_dir, "/", clust_note, "_clust_map_v", as.character(dim(use_dat)[[2]]), "_c", nclust, ".png" )
png(outname, width = 1200, height = 800, res = 200)  # open PNG device
clust_map
dev.off()    


#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done. 
#rmarkdown::render( "Classification_LSSM_PDF.Rmd",   

outname <- paste0( "LSSM_Results_v", as.character(dim(use_dat)[[2]]), "_c", nclust )
rmarkdown::render( "Classification_LSSM_results.Rmd",   
                   output_format = 'pdf_document',
                   output_dir    = results_dir,
                   output_file = outname )

#----
# FIN.






