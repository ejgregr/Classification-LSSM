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
# 2025/07/15: Start work to complete with original Biannucci data and learnings from DFO version.
# 2025/07/21: Completed processing of original Biannucci data to points
#   --> Start major revision to classification, replacing Spectral rasters with DFO point data.


# TO DO: 
#  1. Reproduce (ish) Romina - use as prototype
#  2. Combine attributes and develop update plan
#  4. Spend a TINY bit of time to see if MSEA data can "improve" fit.
#  5. Can repeated clusterings of same data be attached to same RMD report,
#     OR easier by hand, post-hoc?
#######################################################################################

print('Starting Classification  ... LSSM Version 2: Points not rasters')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "classification_functions.R" )

today <- format(Sys.Date(), "%Y-%m-%d")

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
source_dir <- 'C:/Data/SpaceData/Broughton/netcdf'
data_dir   <- 'C:/Data/Git/Classification-LSSM/Data'
results_dir<- 'C:/Data/Git/Classification-LSSM/Results' 

#NOTEs: 
# -For loading of TIFs, see the DFO version of the classification code.
# -Jul 21 version clipped to KD extents

load( file = paste0( source_dir, '/FVCOM_point_data_2025-07-21.rData' ))
str(kd_pts)
fv_dat <- kd_pts

# Now shorten the names ... 
names(fv_dat) <- gsub("salinity", "salt", names(fv_dat))
names(fv_dat) <- gsub("bottom",   "bott", names(fv_dat))
names(fv_dat) <- gsub("surface",  "surf", names(fv_dat))
colnames(fv_dat)

# remove (X,Y) for the moment
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
fv_dat <- fv_dat[ , !grepl("temp_JUL_surf_min", names(fv_dat)) ] # correlated with max
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

dim(fv_dat)

# Correlation analysis done on untransformed predictors. 
# Now work on Normality 


print( "Transform the data  ... ")
PlotHistos( fv_dat, "Untransformed FVCOM predictors" )

# Show predictors with abs(skew) > 1.
skews <- apply(fv_dat, 2, skewness, na.rm = TRUE)
skews[skews > 1]

# MakeMoreNormal() is hard-coded based on trial and error. 
# Adapted to FVCOM data
tfv_dat <- MakeMoreNormal( fv_dat )

PlotHistos( tfv_dat, "Transformed FVCOM predictors" )

# Now scale and center the transformed data.
print("Centering and scaling  ... ")
x <- scale( tfv_dat, center = T,  scale = T )
done_dat <- as.data.frame( x )
PlotHistos( done_dat, "Transformed and scaled FVCOM predictors" )


#---- Final data selection ---- 
# Any sequential dropping of predictors based on analytic results done here
# So all clustering reflects the changes. 

dim( done_dat )
colnames( done_dat )


#---- Part 2 of 3: Cluster number selection ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

#---- Part 2a: Explore number of clusters using Within-sum-of-squares scree plot ----
# Runs kmeans with increasing number of clusters


nclust   <- 18 # number of clusters for scree plotnsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( done_dat, nclust, randomz, imax, 0 )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 7 # the number of clusters based on scree plot, above.
#nsample <- 500000 # a larger sample for more robust classification

#sidx <- sample( 1:length( stack_data_clean[ , 1] ), nsample )
#samp <- stack_data_clean[ sidx, ]

samp <- done_dat

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

cs <- cluster_result$cluster
c_dist <- dist(samp)
sk <- silhouette(cs, c_dist)

plot(sk, col = 1:nclust, border=NA, main = "" )

#---- Part 3: Detailed examination of N clusters  ----
#---- Part 3a: Show cluster groupings using PCA ----

z <- as.data.frame( cbind(cluster = cs, done_dat ) )
pca_results <- ClusterPCAall( z )
names( pca_results ) <- c("loadings", "plot1","plot2")

pca_results$plot1
#Percentage of variance explained by dimensions
#eigenvalue <- round(get_eigenvalue(res_pca), 1)
#var_percent <- eigenvalue$variance.percent

#---- Part 3b: Violins of predictor contributions to WORKING clusters ----

x <- done_dat
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
# Jul22: Currently using points here. 
#   For a raster, see the code in the DFO classification version.

pts <- data.frame(x = kd_pts$X, y = kd_pts$Y, cluster = cluster_result$cluster)

# 32609 is the CRS for UTM Zone 9N, the projection of the FVCOM data
sf_pts <- st_as_sf(pts, coords = c("x", "y"), crs = 32609)

hist( pts$cluster )
ggplot(sf_pts, aes(color = factor(cluster))) + geom_sf()

#Output a shape file
st_write(sf_pts, "cluster_points.shp", delete_layer = TRUE)




# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent") # Max for Accent is 8



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






