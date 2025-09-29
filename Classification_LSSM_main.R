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
# 2025/09/10: Developed a seasonal plan with Chris. Optimistic this will lead to a conclusion.

# To DO:
# - Examine predictors with cor = 1.0. FVCOM or data transfer problem?

################################################################################

print('Starting Classification  ... LSSM Version 2: Points not rasters')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
 source( "classification_functions.R" )

# Sourcing here reloads FVCOM data into kd_pts. Includes reading NETCDFs, interpolations, and clipping.
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
# Sept 5 version. Minor code tweaks
# Loads structure kd_pts, 
load( file = paste0( data_dir, '/FVCOM_point_data_2025-09-11.rData' ))

str(kd_pts)
fv_dat <- kd_pts

# Shorten some names ... 
names(fv_dat) <- gsub("salinity", "salt", names(fv_dat))
names(fv_dat) <- gsub("bottom",   "bott", names(fv_dat))
names(fv_dat) <- gsub("surface",  "surf", names(fv_dat))

# And re-order the column names to facilitate comparison of figures
# Nice bit of perl from chatgpt! 
colnames(fv_dat) <- sub(
      "(^[^_]+)_([A-Z]+)_(.*)",
      "\\L\\2\\E_\\U\\1\\E_\\3",
      colnames(fv_dat),
    perl = TRUE
  )


# remove (X,Y) for the analysis
#plot( kd_pts$X, kd_pts$Y)
fv_dat    <- fv_dat[ , !(names(fv_dat) %in% c("X", "Y" )) ]
fv_coords <- kd_pts[ , c("X", "Y" ) ]

#---- Correlation and standardization section ----

### Reduce correlation for July predictors
jul_dat <- fv_dat[ , !grepl("dec", names(fv_dat)) ]

# Pass 1: remove all cors > 0.9
showHighestPairs( jul_dat, 0.9)
jul_dat <- jul_dat[ , !grepl("jul_SALT_bott_ave", names(jul_dat)) ]
jul_dat <- jul_dat[ , !grepl("jul_TEMP_bott_ave", names(jul_dat)) ]
jul_dat <- jul_dat[ , !grepl("jul_SALT_surf_ave", names(jul_dat)) ]
jul_dat <- jul_dat[ , !grepl("jul_TEMP_surf_ave", names(jul_dat)) ]
jul_dat <- jul_dat[ , !grepl("jul_TEMP_surf_min", names(jul_dat)) ]

# Pass 2: remove all cors > 0.8
showHighestPairs( jul_dat, 0.8)
# Remove all July bottom currents
jul_dat <- jul_dat[ , !grepl("jul_CS_bott", names(jul_dat)) ]

jul_dat <- jul_dat[ , !grepl("jul_CS_surf_min", names(jul_dat)) ]
jul_dat <- jul_dat[ , !grepl("jul_SALT_bott_min", names(jul_dat)) ]

# Pass 3: remove all cors > 0.5
# NOTE cor=0.7 and 0.6 fail cuz there is only 1 ... 
showHighestPairs( jul_dat, 0.4)
# Cor btwn max bott temp and max bottom salt = -0.72 but both retain cuz chemistry.
names(jul_dat)


### Reduce correlation for December predictors
dec_dat <- fv_dat[ , !grepl("jul", names(fv_dat)) ]

# Pass 1: remove all cors > 0.9
showHighestPairs( dec_dat, 0.9)
# reduce chemistry correlations
dec_dat <- dec_dat[ , !grepl("dec_SALT_bott_ave", names(dec_dat)) ]
dec_dat <- dec_dat[ , !grepl("dec_TEMP_bott_ave", names(dec_dat)) ]
dec_dat <- dec_dat[ , !grepl("dec_SALT_surf_ave", names(dec_dat)) ]
dec_dat <- dec_dat[ , !grepl("dec_TEMP_surf_ave", names(dec_dat)) ]

# Pass 2: remove all cors > 0.8
showHighestPairs( dec_dat, 0.8)
# Remove all Dec bottom currents
dec_dat <- dec_dat[ , !grepl("dec_CS_bott", names(dec_dat)) ]
dec_dat <- dec_dat[ , !grepl("dec_CS_surf_min", names(dec_dat)) ]
# Remove min of temp correlation
dec_dat <- dec_dat[ , !grepl("dec_TEMP_bott_min", names(dec_dat)) ]

# Pass 3: remove all cors > 0.7
showHighestPairs( dec_dat, 0.7)
# Remove min of salt ... 
dec_dat <- dec_dat[ , !grepl("dec_SALT_bott_min", names(dec_dat)) ]
# Remove max of temp... 
dec_dat <- dec_dat[ , !grepl("dec_TEMP_surf_max", names(dec_dat)) ]

# Pass 3: remove all cors > 0.6
showHighestPairs( dec_dat, 0.4)

names(dec_dat)

# Prepare histograms of seasonal Untransformed, uncorrelated FVCOM predictors
hist1 <- PlotHistos( jul_dat, "Original FVCOM predictor values" )
hist2 <- PlotHistos( dec_dat, "Original FVCOM predictor values" )

#---- Transformation and scaling of uncorrelated predictors  ----
print( "Transform the data  ... ")

# As above, we'll do summer then winter ... 

# JUL predictors with abs(skew) > 1.
jul_skew <- apply(jul_dat, 2, skewness, na.rm = TRUE)
jul_skew[abs(jul_skew) > 1]

# DEC predictors with abs(skew) > 1.
dec_skew <- apply(dec_dat, 2, skewness, na.rm = TRUE)
dec_skew[abs(dec_skew) > 1]

# Only one predictor from each season exceeds a skew of 2.0
# Those are fixed in place here. MakeMoreNormal() not required.
x <- jul_dat[, "jul_SALT_surf_min"]
floor <- 21
y <- ifelse(x < floor, floor, x)
y <- y^3
jul_dat[, "jul_SALT_surf_min"] <- y

x <- dec_dat[, "dec_SALT_surf_min"]
floor <- 21
y <- ifelse(x < floor, floor, x)
y <- y^8
dec_dat[, "dec_SALT_surf_min"] <- y

#skewness(y)
#hist(y, breaks=25)

# Now scale and center the transformed data.
print("Centering and scaling  ... ")

x <- scale( jul_dat, center = T,  scale = T )
jul_dat <- as.data.frame( x )

x <- scale( dec_dat, center = T,  scale = T )
dec_dat <- as.data.frame( x )

hist3 <- PlotHistos( jul_dat, "Scaled and transformed predictors" )
hist4 <- PlotHistos( dec_dat, "Scaled and transformed predictors" )


#---- Final data selection ---- 
# Any sequential dropping of predictors based on analytic results done here
# 09/12: Not sure we be doing any of this now with seasonal comparison. 


#---- Cluster analysis Part 1 - Number of clusters  ----
  set.seed <- 42 # Seed for reproducibility
  randomz  <- 20 # the number of randomizations for kmeans to do.
  imax     <- 25 # maximum iterations to try for convergence
  
  # Explore number of clusters using Within-sum-of-squares scree plot
  # Runs kmeans with increasing number of clusters
  nclust     <- 18 # number of clusters for scree plotnsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
  scree_plot1 <- MakeScreePlot( jul_dat, nclust, randomz, imax, 0 )
  scree_plot2 <- MakeScreePlot( dec_dat, nclust, randomz, imax, 0 )
  
  
  #---- Cluster analysis Part 2 - Create working set of N clusters ----
  
  #-- SETUP for clustering and figure labeling ----
  nclust  <- 9 # the number of clusters based on scree plot, above.
  cl_text <- "9"
  
  #-- RUn the clusters
  names( jul_dat )
  names( dec_dat )
  
  jul_clust <- kmeans(jul_dat, centers=nclust, nstart=randomz, iter.max=imax) 
  dec_clust <- kmeans(dec_dat, centers=nclust, nstart=randomz, iter.max=imax) 
  

#---- Part 3: Cluster diagnostic plots  ----

# The cluster and data to use in the plots
#---- Heat map of within-cluster standard deviations ----
heat_map1 <- PlotHeatMap( jul_clust$cluster, jul_dat ) 
heat_map2 <- PlotHeatMap( dec_clust$cluster, dec_dat ) 

#---- Examine silhouette plots  ----
# This code duplicated in RMD file so plots are generated. 

c_dist <- dist( jul_dat )
sk <- silhouette(jul_clust$cluster, c_dist)
# Call generic plot() in the r script
jul_silplot <- plot(sk, col = 1:nclust, border=NA, main = "" )

c_dist <- dist( dec_dat )
sk <- silhouette(dec_clust$cluster, c_dist)
dec_silplot <- plot(sk, col = 1:nclust, border=NA, main = "" )


#---- Show cluster groupings using PCA ----
jul_pca <- ClusterPCAall( cbind(cluster = jul_clust$cluster, jul_dat ) )
names( jul_pca ) <- c("loadings", "plot1","plot2")

dec_pca <- ClusterPCAall( cbind(cluster = dec_clust$cluster, dec_dat ) )
names( dec_pca ) <- c("loadings", "plot1","plot2")


#Percentage of variance explained by dimensions
jul_eigens <- round(get_eigenvalue(jul_pca[[1]]), 2)
dec_eigens <- round(get_eigenvalue(dec_pca[[1]]), 2)

var_percent <- jul_eigens$variance.percent

#---- Violins of predictor contributions to WORKING clusters ----
jul_viols <- PlotViolins( jul_clust$cluster, jul_dat )
dec_viols <- PlotViolins( dec_clust$cluster, dec_dat )


#---- Part 4: Map the working point clusters ----
# NB: Using points here. For a raster, see the DFO version.

jul_map <- PlotMap( fv_coords, jul_clust, "July Clusters" )
dec_map <- PlotMap( fv_coords, dec_clust, "December Clusters" )

jul_map / dec_map

# NEXT: Standardize the cluster numbers by space


#================ STOP HERE for now ====================


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

outname <- paste0( "LSSM_Results_v", as.character(dim(dec_dat)[[2]]), "_c", nclust )
rmarkdown::render( "Classification_LSSM_results.Rmd",   
                   output_format = 'pdf_document',
                   output_dir    = results_dir,
                   output_file = outname )

#----
# FIN.






