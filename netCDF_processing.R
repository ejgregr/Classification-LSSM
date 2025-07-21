#----------------------------------------------------------------------------
# Script:  NetCDF_processing.R
# Created: July 2025 with help from ChatGPT to process files provided by Laura Bianucci.
#
# Notes:
#  - 2025/07/16: 
#
#================================== Load require packages ======================

# check for any required packages that aren't installed and install them
required.packages <- c( "ggplot2", "ncdf4", "terra", "akima", "sf", "R.matlab", "clue", "viridis")

uninstalled.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]

# install any packages that are required and not currently installed
if(length(uninstalled.packages)) install.packages(uninstalled.packages)

# require all necessary packages
lapply(required.packages, require, character.only = TRUE)
#lapply(required.packages, library, character.only = TRUE)


#=========================== Data sources and constants =====================================

# proj4 string : Albers projection, NAD83 datum
albers_crs <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

nc_dir      <- 'C:/Data/SpaceData/Broughton/netcdf' 
results_dir <- 'C:/Data/Git/Classification-LSSM/Results' 

#============================== Functions ======================================

# Pull all the variables from a netcdf file. 
GetNCVariables <- function(ncfile) {
  # Open NetCDF file
  nc <- nc_open(ncfile)
  
  # Get all variable names (excluding dimension variables)
  varnames <- names(nc$var)
  
  # Create named list of variables
  vars <- lapply(varnames, function(v) ncvar_get(nc, v))
  names(vars) <- varnames
  
  # Close file
  nc_close(nc)
  
  # Return list of variables
  return(vars)
}

# Parse filenames into relevant bits. 
ParseNCFilename <- function(filename) {
  # Example: QCS_DEC2018N_1month_TSuv_bottom_ave.nc
  month <- sub(".*_(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)[0-9]{4}.*", "\\1", filename)
  level <- ifelse(grepl("bottom", filename, ignore.case = TRUE), "bottom", "surface")
  metric <- ifelse(grepl("min", filename), "min",
                   ifelse(grepl("max", filename), "max", "ave"))
  list(month = month, level = level, metric = metric)
}

# Pull the desired variables from a specified file
ExtractVars <- function(ncfile) {
  nc <- nc_open(ncfile)
  vars <- list(
    salinity = as.vector(ncvar_get(nc, "salinity")),
    temp     = as.vector(ncvar_get(nc, "temp")),
    u        = as.vector(ncvar_get(nc, "u")),
    v        = as.vector(ncvar_get(nc, "v")),
    x        = as.vector(ncvar_get(nc, "x")),
    y        = as.vector(ncvar_get(nc, "y"))
  )
  nc_close(nc)
  return(vars)
}

# Interpolate u and v values from faces to nodes using 
# n2f: list where node2faces[[i]] = indices of elements for node i
# uv_data: named list or data frame of element-based vectors
# --> all same length, names=column names)
InterpElemToNode <- function(n2f, uv_data) {
  uv_data <- as.list(uv_data)  # Accept data frame or list
  nnode   <- length(n2f)
  res     <- list()
  
  for (varname in names(uv_data)) {
    values <- uv_data[[varname]]
    node_values <- numeric(nnode)
    for (i in seq_len(nnode)) {
      idx <- n2f[[i]]
      node_values[i] <- mean(values[idx], na.rm = TRUE)
    }
    res[[varname]] <- node_values
#    print(paste("Processing variable:", varname))
#    print(summary(values))
  }
  return( as.data.frame(res))
}


# Return current speed for each pair of u,v (of the form 'u_DEC_bottom_mean').
# Expecting 12 columns from an input of 24.  
CalcCurspd <- function(df) {
  # Find all u columns
  u_cols <- grep("^u_", names(df), value = TRUE)
  
  # For each u column, find the corresponding v column
  curspd_list <- list()
  for (u_col in u_cols) {
    # Find corresponding v column by replacing u_ with v_
    v_col <- sub("^u_", "v_", u_col)
    if (v_col %in% names(df)) {
      # Compute current speed
      curspd <- sqrt(df[[u_col]]^2 + df[[v_col]]^2)
      # Build new column name
      curspd_col <- sub("^u_", "CS_", u_col)
      curspd_list[[curspd_col]] <- curspd
    }
  }
  # Return as data frame with columns for each speed
  as.data.frame(curspd_list)
}


# A kludge to temporarily deal with the node/faces problem
TruncToShortest <- function(lst) {
  minlen <- min(lengths(lst))
  lapply(lst, function(x) x[seq_len(minlen)])
}


# Interpolates a vector of values at x/y nodes to a raster grid
# Keeping just for fun. 
InterpToGrid <- function(x, y, values, nx, ny) {
  ok <- !(is.na(x) | is.na(y) | is.na(values))
  interp_res <- akima::interp(x[ok], y[ok], values[ok], nx = nx, ny = ny, duplicate = "mean")
  rast(t(interp_res$z), extent = c(range(interp_res$x), range(interp_res$y)))
}

#============================== Main Code ======================================

#---- Load predictors from nc files ----
# make a list of predictors from FVCOM data 
nc_list <- list.files(path = nc_dir, pattern = '\\.nc$', full.names = TRUE)

# Examining the netCDF contents  ... 
vars <- GetNCVariables(  nc_list[1] )
str(vars, max.level=1)

# Each file contains salinity, temp, u, and v. 
# Month (JUL, DEC), level (bottom, surface), and metric (min, max, ave) are contained in the files names. 

# Initialize target structure
all_data <- list()

for (file in nc_list) {
  attrs <- ParseNCFilename(file)
  varlist <- ExtractVars(file) # Each variable is a vector (length = n nodes)
  
  for (varname in names(varlist)) {
    colname <- paste(varname, attrs$month, attrs$level, attrs$metric, sep = "_")
    all_data[[colname]] <- varlist[[varname]]
  }
}
# NOTE: Currents (u,v) have a different length because they are on faces. 
#--- DATA:
# 123,598 appears to be the number of nodes
# 221,232 appears to be the number of faces

#---- Interpolate currents to nodes ----

# Start with the connectivity matrix from Laura's .mat file. 
matfile <- readMat( paste0( nc_dir, "/QCS_GridInfo.mat") )
c_table <- matfile$nv
dim(c_table) 

# Now build a simple mean of all the u,v values on the faces that connect to each node. 
nnode <- max(c_table)
nele  <- nrow(c_table)

# Initialize node to faces list
node2faces <- vector("list", nnode)

# build a node to faces list
for (elem in 1:nele) {
  for (j in 1:3) {
    node <- c_table[elem, j]
    node2faces[[node]] <- c(node2faces[[node]], elem)
  }
}

# Now do the interplation 
uv_names <- grep("^[uv]", names(all_data), value = TRUE, ignore.case = TRUE)
uv_data  <- all_data[ uv_names ]
uv_interp <- InterpElemToNode( node2faces, uv_data )
names(uv_interp)

# And finally combine u and v into a current speed
uv_cspd <- CalcCurspd( uv_interp )
#head(uv_cspd)

# Now combine the original node-based values with the calculated current speed 
other_names <- setdiff( names(all_data), uv_names )
other_data  <- all_data[ other_names ]

big_df <- as.data.frame( c( other_data, uv_cspd ))

# Combine all columns into a data frame ... this would be nice but ... 
# measures are on nodes and u,v are on the faces. 
str(big_df)

### Now have 60 variables from the NC files all the same length (123,598)

#---- Final data preparation step ----
# Refactor the dataframe to remove redundancies
# First pull the coordinates - only need them once! :)
xy <- data.frame( x = big_df$x_DEC_bottom_ave, y = big_df$y_DEC_bottom_ave )
# Now drop all the coord and current fields - currents will be handled separately
big_df <- big_df[ , !grepl("^(x|y)", names(big_df))]
# Reattach the coordinates
big_df <- cbind( xy, big_df )
dim(big_df)

names(big_df)
#str(big_df)

#---- Clip to Klokwade Damsxi ----
# We're going to work in UTM Zone 9 so that the clipped points can be used 
# to then trim the rest of the data. 

kd_bounds <- st_read( "c:\\Data\\SpaceData\\Broughton\\FVCOM_trim.shp")
utm_pts <- st_as_sf(big_df, coords = c("x", "y"), crs = "+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
albers_pts <- st_transform(utm_pts, coords = c("x", "y"), crs = albers_crs)
kd_pts <- albers_pts[ kd_bounds, ]

ggplot() +
  geom_sf(data = kd_bounds, fill = NA, color = "blue", lwd = 1) +
  geom_sf(data = kd_pts, color = "red", size = 1) +
  theme_minimal()

# turn kd_pts back into a df after messing about with sf points for the clip.
coords <- st_coordinates(kd_pts)
kd_pts <- cbind( coords,st_drop_geometry(kd_pts) )
str(kd_pts)


# Confirm u/v processing look ok ... 
ggplot(kd_pts, aes(x = X, y = Y, color = CS_DEC_bottom_ave)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()


#-------------- Correlation bits ------------------
# Copy cor code from main script ... 
#---- Correlation across data layers 
cor_table <- round( cor( sum_dat ), 3)
cor_table[lower.tri(cor_table, diag=TRUE)] <- NA

high_rows <- apply(cor_table, 1, function(row) any(row > 0.6, na.rm = TRUE))
z <- cor_table[ high_rows, ]
z
# Start removing variables. 
sum_dat <- sum_dat[ , !grepl("salinity_JUL_bottom_ave", names(sum_dat))]
sum_dat <- sum_dat[ , !grepl("temp_JUL_bottom_ave", names(sum_dat))]
sum_dat <- sum_dat[ , !grepl("salinity_JUL_surface_ave", names(sum_dat))]
sum_dat <- sum_dat[ , !grepl("temp_JUL_surface_ave", names(sum_dat))]
sum_dat <- sum_dat[ , !grepl("temp_JUL_surface_min", names(sum_dat))]
sum_dat <- sum_dat[ , !grepl("salinity_JUL_surface_max", names(sum_dat))]

# This is the only one proving a bit difficult cuz both min and max are correlated
# with salinity_JUL_surface_min. Maybe leave with this correlation for now.

#sum_dat <- sum_dat[ , !grepl("salinity_JUL_bottom_min", names(sum_dat))]
dim(sum_dat)


#---- Initial clustering ----

### WORKING WITH ALL PREDICTORS - UNTIL CAN INCLUDE CURRENTS FOR COR WORK ###

# Lets look at seasonal differences ... 
sum_dat <- kd_pts[ ,grepl( "(X|Y|JUL)", names(kd_pts) )]
win_dat <- kd_pts[ ,grepl( "(X|Y|DEC)", names(kd_pts) )]

set.seed(43)
k <- 6

# Make clusters but without (X, Y)
sum_clust <- kmeans( sum_dat[ , !(names(sum_dat) %in% c("X", "Y")) ], centers = k )
win_clust <- kmeans( win_dat[ , !(names(win_dat) %in% c("X", "Y")) ], centers = k )

sum_dat$cluster <- sum_clust$cluster
win_dat$cluster <- win_clust$cluster

# Aligning cluster labels ... 
clustperm <- clue::solve_LSAP( table(sum_clust$cluster, win_clust$cluster), maximum = TRUE)
# Apply permutation to clusteringB to best match clusteringA
win_dat$cluster <- clustperm[ win_clust$cluster ]

#---- Plotting
library(ggplot2)
a <- ggplot(sum_dat, aes(x = X, y = Y, color = factor(cluster))) +
  geom_point(size = 1) +
  coord_equal() +
  theme_minimal()
a
ggsave( paste0( results_dir, "/summer_6cluster.png"), plot = a, width = 6, height = 4, dpi = 300)


library(ggplot2)
a <- ggplot(win_dat, aes(x = X, y = Y, color = factor(cluster))) +
  geom_point(size = 1) +
  coord_equal() +
  theme_minimal()
a
ggsave( paste0( results_dir, "/winter_6cluster.png"), plot = a, width = 6, height = 4, dpi = 300)





#---------------------- JOINT CLUSTERING?? -------------------
names(sum_dat) <- gsub("_JUL_", "_", names(sum_dat))
names(win_dat) <- gsub("_DEC_", "_", names(win_dat))
sum_dat$season <- "summer"
win_dat$season <- "winter"

combined <- rbind(sum_dat, win_dat)
# Only cluster on environmental columns
env_cols <- setdiff(names(combined), c("X", "Y", "season"))

k <- 6
set.seed(43)
joint_clust <- kmeans(combined[ , env_cols ], centers = k)

combined$cluster <- joint_clust$cluster

# Split back out
sum_dat$cluster <- combined$cluster[combined$season == "summer"]
win_dat$cluster <- combined$cluster[combined$season == "winter"]

#---- Plotting
a <- ggplot(sum_dat, aes(x = X, y = Y, color = factor(cluster))) +
  geom_point(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +  # controls legend point size
  coord_equal() +
  theme_minimal()
a
ggsave( paste0( results_dir, "/summer_6cluster.png"), plot = a, width = 6, height = 4, dpi = 300)


a <- ggplot(win_dat, aes(x = X, y = Y, color = factor(cluster))) +
  geom_point(size = 1) +
  coord_equal() +
  theme_minimal()
a
ggsave( paste0( results_dir, "/winter_6cluster.png"), plot = a, width = 6, height = 4, dpi = 300)













#---------------------- Vestigial code -------------------------------------
#  Interpolation, projection, and rasterization of netCDF data

# Example usage for FVCOM node data:
# Suppose you want surface temperature at time 1, surface layer (siglay = 1):

x <- vars$x   # node x-coordinates
y <- vars$y   # node y-coordinates
xres = 200    # raster resolution
yres = 200    # raster resolution
temp <- vars$temp  # [node, siglay, time] array


# Interpolate to regular grid and create raster
r_temp <- InterpToGrid( vars$x, vars$y, vars$h, xres, yres )

# Flip to correct orientation on y axis - apparetly a common problem with 
# netcdf to raster conversion.
r_temp <- flip(r_temp, direction = "vertical")


# Assign coordinate system if you know it (example is UTM zone 9N)
crs(r_temp) <- "+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs"
r_temp <- project( r_temp, spat_ref )
r_temp <- raster::crop( r_temp, amask )

plot(r_temp, main = "Monthly mean temperature (FVCOM)")

writeRaster(r_temp, paste0( data_dir, "/fvcom_temp2.tif"), overwrite = TRUE)


# Assemble the coordinates for display in Arc 
# Goal was to try Python-based barrier-aware interpolation until found the 
# lack of data in KKD. 
x <- cbind( vars$x, vars$y, vars$h )
colnames(x) <- c( "x", "y", "h")
head(x)
write.csv( x, file = paste0( data_dir, '/FVCOM_depth_points.txt' ), row.names = FALSE)
