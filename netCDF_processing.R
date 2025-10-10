#----------------------------------------------------------------------------
# Script:  NetCDF_processing.R
# Created: July 2025 with help from ChatGPT to process FVCOM files provided by Laura Bianucci.
#
# Notes:
#  - 2025/09/09: Largely complete and stand-alone. 
#     Handles NETCDF loading, projecting, and clipping to spatial extents defined using a shape file.
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


#=============== Directories Will be created if they don't exist. ==============
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
vars <- GetNCVariables( nc_list[1] )
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

# Now do the interpolation 
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
# Re-factor the dataframe to remove redundancies
# First pull the coordinates - only need them once! :)
xy <- data.frame( x = big_df$x_DEC_bottom_ave, y = big_df$y_DEC_bottom_ave )
# Now drop all the coord and current fields - currents will be handled separately
big_df <- big_df[ , !grepl("^(x|y)", names(big_df))]
# Reattach the coordinates
big_df <- cbind( xy, big_df )
dim(big_df)

names(big_df)
str(big_df)

#---- Clip to specified spatial extents (kd_bounds) ----
# We're going to work in UTM Zone 9 so that the clipped points can be used 
# to then trim the rest of the data. 

# Assign the NAD83 UTM projection
utm_pts <- st_as_sf(big_df, coords = c("x", "y"), crs = UTM_crs )
# Transform to BC Albers to match rest of data
albers_pts <- st_transform(utm_pts, albers_crs)

kd_pts <- albers_pts[ kd_bounds, ]
coast <- st_intersection( coast, kd_bounds )
# OR just use the unclipped data ... 
# kd_pts <- albers_pts

ggplot() +
  geom_sf(data = kd_bounds, fill = NA, color = "blue", lwd = 1) +
  geom_sf(data = kd_pts, color = "red", size = 1) +
  geom_sf(data = coast, fill = NA, color = "black", lwd = 0.8) +
  theme_minimal()

# turn kd_pts back into a df after messing about with sf points for the clip.
coords <- st_coordinates(kd_pts)
kd_pts <- cbind( coords, st_drop_geometry(kd_pts) )
str(kd_pts)


# Confirm u/v processing look ok ... 
ggplot(kd_pts, aes(x = X, y = Y, color = CS_DEC_bottom_ave)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()


# Save the results!
target_dir <- data_dir # From main script.
save( kd_pts, file = paste0( target_dir, '/FVCOM_point_data_', today, '.rData' ))


### FIN!
