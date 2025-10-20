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
# 2025/10/01: Updated data with Laura and predictors with cor = 1.0 gone.
#   --> Settled on expanded area to include the entire 'estuary' in the classification.
#   --> Implement approach to have extent-specific predictor prep.
#   --> Good progress on RMD file.
#   --> Have prepped predictors for Kwakwa region
#
# TO DO:- Standardize palette across clusters in results figs
#       - Add coastline to map figures, and make them as large as possible
#       - Fix Jul transforms for KwaKwa
#       - Add a) and b) to relevant figure captions
#       - Compare Kwakwa Jul PCA with remaining correlations in the data
################################################################################

print('Starting Classification  ... LSSM Version 2: Points not rasters')
rm(list = ls(all = T))  # Erase environment.
today <- format(Sys.Date(), "%Y-%m-%d")

#======================== Directories and constants ===========================
# Projections as EPSG codes
albers_crs <- 3005 # Or for newer datasets: albers_crs <- 3153
UTM_crs    <- 26909 # For Zone 9N NAD83. Or for WGS84: 32609

# Will be created if they don't exist.
source_dir <- 'C:/Data/SpaceData/Broughton/netcdf'
data_dir   <- 'C:/Data/Git/Classification-LSSM/Data'
results_dir <- 'C:/Data/Git/Classification-LSSM/Results'

# Load necessary packages and functions ...
source("classification_functions.R")

# Sourcing here reloads FVCOM data into kd_pts. Includes reading NETCDFs, interpolations, and clipping.
#source("netCDF_processing.R")

# NOTES:
# - For loading of TIFs, and the necessary subsampling, see the DFO version of the classification code.
# - The FVCOM data are only sub-sampled spatially using a shapefile from ArcGIS.

#---- Load FVCOM point data ----
# Loads kd_pts dataframe
load(file = paste0(data_dir, '/FVCOM_point_data_2025-10-09_kwakwa.rData'))
#load(file = paste0(data_dir, '/FVCOM_point_data_2025-09-11_village.rData'))

# Load shape files for FCVOM bounding box (for clipping) and coastline (for plotting)
kd_bounds <- st_read( "c:\\Data\\SpaceData\\Broughton\\FVCOM_trim.shp")
coast     <- st_read( "c:\\Data\\SpaceData\\Broughton\\coastline.shp")

# cLIP coast to bounds
coast <- st_intersection( coast, kd_bounds )

str(kd_pts)
fv_dat <- kd_pts

# Shorten some names ...
names(fv_dat) <- gsub("salinity", "salt", names(fv_dat))
names(fv_dat) <- gsub("bottom", "bott", names(fv_dat))
names(fv_dat) <- gsub("surface", "surf", names(fv_dat))

# And re-order the column names to facilitate comparison of figures
# Nice bit of perl from chatgpt!
colnames(fv_dat) <- sub("(^[^_]+)_([A-Z]+)_(.*)",
                        "\\L\\2\\E_\\U\\1\\E_\\3",
                        colnames(fv_dat),
                        perl = TRUE)

# remove (X,Y) for the analysis
#plot( kd_pts$x, kd_pts$Y)
fv_dat    <- fv_dat[, !(names(fv_dat) %in% c("X", "Y"))]
fv_coords <- kd_pts[, c("X", "Y")]


#---- Correlation and standardization section ----
# Specific to extents of FCVOM used.
# FUNCTIONS exist for each extents to preserve the parameter selection
# and transforms applied. These include:
# ParamVilliageSea  - The centre of Broughton. Incl. Village sea
# ParamKwakwakuitil - FVCOM subset for full Broughton estuary
# ParamFVCOM        - The full FVCOM extents.

# THE TOP DOWN approach focuses on removing correlations, with priority variables. 
#source( 'ParamVillageSeaJul.R' )
#source( 'ParamVillageSeaDec.R' )
#source( 'ParamKwakwakuitilJul_topdown.R' )
#source( 'ParamKwakwakuitilDec_topdown.R' )

# BOTTOM UP approach starts with what we expect might be ecologically reasonable.


#---- Process July predictors ----
jul_dat <- fv_dat[, !grepl("dec", names(fv_dat))]
source( 'ParamKwakwakuitilJul_botup.R' )

#---- Process December predictors ----
dec_dat <- fv_dat[, !grepl("jul", names(fv_dat))]
source( 'ParamKwakwakuitilDec_botup.R' )

hist1 / hist3
hist2 / hist4

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
scree_plot1 <- MakeScreePlot(jul_dat, nclust, randomz, imax, 0)
scree_plot2 <- MakeScreePlot(dec_dat, nclust, randomz, imax, 0)


#---- Cluster analysis Part 2 - Create working set of N clusters ----

#-- SETUP for clustering and figure labeling ----
nclust  <- 4 # the number of clusters based on scree plot, above.
cl_text <- "4"

#-- RUn the clusters
names(jul_dat)
names(dec_dat)

jul_clust <- kmeans(jul_dat,
                    centers = nclust,
                    nstart = randomz,
                    iter.max = imax)
dec_clust <- kmeans(dec_dat,
                    centers = nclust,
                    nstart = randomz,
                    iter.max = imax)

#---- Standardize pairs of clusters ----
# ---- Usage ----
# Assume you already have:
#   km_ref   <- kmeans(X_ref, centers = K, nstart = 25)
#   km_new   <- kmeans(X_new, centers = K, nstart = 25)
# and a coords data frame `pts` with columns x,y in the SAME row order as X_ref/X_new.

km_ref <- jul_clust
km_new <- dec_clust
align <- AlignClusters(km_ref, km_new)

# Cross-tab BEFORE ahd after alignment
#cat("\nContingency (ref vs new) BEFORE alignment:\n")
table(km_ref$cluster, km_new$cluster)
table(km_ref$cluster, align$aligned)

cat("\nMapping (index = new label, value = ref label):\n")
print(align$mapping)

cat("\nContingency used for assignment:\n")
print(align$contingency)

# Check AFTER alignment
cat("\nContingency (ref vs aligned) AFTER alignment:\n")
print(table(km_ref$cluster, align$aligned))

# jul_clust is unchanged.
dec_clust <- ApplyAlignment(dec_clust, align$mapping)


# Create a consistent palette with vivid, distinct colors (HCL space)
random_pal <- function(K, seed = NULL, l = 65, c = 100) {
  if (!is.null(seed))
    set.seed(seed)
  hues <- sample(seq(0, 360 - 360 / K, length.out = K))
  stats::setNames(grDevices::hcl(h = hues, c = c, l = l), as.character(1:K))
}

our_pal <- random_pal( nclust, seed=13 )

#---- Part 3: Cluster diagnostic plots  ----

# The cluster and data to use in the plots
#---- Heat map of within-cluster standard deviations ----
heat_map1 <- PlotHeatMap(jul_clust$cluster, jul_dat)
heat_map2 <- PlotHeatMap(dec_clust$cluster, dec_dat)

#---- Examine silhouette plots  ----
# This code duplicated in RMD file so plots are generated.

c_dist <- dist(jul_dat)
sk <- silhouette(jul_clust$cluster, c_dist)
# Call generic plot() in the r script
plot(sk,
     col = our_pal,
     border = NA,
     main = "")

c_dist <- dist(dec_dat)
sk <- silhouette(dec_clust$cluster, c_dist)
plot(sk,
     col = our_pal,
     border = NA,
     main = "")


#---- Show cluster groupings using PCA ----
jul_pca <- ClusterPCAall(cbind(cluster = jul_clust$cluster, jul_dat))
names(jul_pca) <- c("loadings", "plot1", "plot2")

dec_pca <- ClusterPCAall(cbind(cluster = dec_clust$cluster, dec_dat))
names(dec_pca) <- c("loadings", "plot1", "plot2")


#Percentage of variance explained by dimensions
jul_eigens <- round(get_eigenvalue(jul_pca[[1]]), 2)
dec_eigens <- round(get_eigenvalue(dec_pca[[1]]), 2)

var_percent <- jul_eigens$variance.percent

#---- Violins of predictor contributions to WORKING clusters ----
jul_viols <- PlotViolins(jul_clust$cluster, jul_dat, our_pal)
dec_viols <- PlotViolins(dec_clust$cluster, dec_dat, our_pal)


#---- Part 4: Map the working point clusters ----
# NB: Using points here. For a raster, see the DFO version.

x <- as.factor( jul_clust$cluster )
names( our_pal ) <- levels( x )
theme_set(theme_gray()) # Reset ggplot() defaults

# uses default our_pal
jul_map <- PlotMap(fv_coords, jul_clust$cluster, "July Clusters", our_pal)
dec_map <- PlotMap(fv_coords, dec_clust$cluster, "December Clusters", our_pal)


#---- Knit and render Markdown file to PDF -----
# First had to install the library tinytex.
# then run >tinytex::install_tinytex()
# ... and done.
#rmarkdown::render( "Classification_LSSM_PDF.Rmd",

outname <- paste0("LSSM_Results_v", as.character(dim(dec_dat)[[2]]), "_c", nclust)
rmarkdown::render(
  "Classification_LSSM_results.Rmd",
  output_format = 'pdf_document',
  output_dir    = results_dir,
  output_file = outname
)

#================ STOP HERE for now ====================

#----  Output a shape file ----
# sf_pts created above
outname <- paste0("LSSM_Results_v", as.character(dim(use_dat)[[2]]), "_c", nclust)
st_write(sf_pts, paste0( results_dir, outname, ".shp"), delete_layer = TRUE)

# Define color palette - Max N for Accent is 8
#pal_clust <- brewer.pal(n = nclust, "Accent") 


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


#----
# FIN.



