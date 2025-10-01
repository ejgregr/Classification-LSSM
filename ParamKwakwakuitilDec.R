#---- Section ParamKwakwakuitilDec ----

# Pass 1: remove all cors > 0.9
showHighestPairs(dec_dat, 0.9)
# reduce chemistry correlations
dec_dat <- dec_dat[, !grepl("dec_SALT_bott_max", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_SALT_bott_min", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_TEMP_surf_min", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_TEMP_bott_min", names(dec_dat))]

dec_dat <- dec_dat[, !grepl("dec_SALT_bott_ave", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_SALT_surf_ave", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_TEMP_surf_max", names(dec_dat))]

# Pass 2: remove all cors > 0.8
showHighestPairs(dec_dat, 0.8)
# Remove all Dec bottom currents
dec_dat <- dec_dat[, !grepl("dec_TEMP_bott_max", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_SALT_surf_max", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_CS_bott_min", names(dec_dat))]

showHighestPairs(dec_dat, 0.5)

# Prepare histograms of seasonal Untransformed, uncorrelated FVCOM predictors
hist2 <- PlotHistos(dec_dat, "Original FVCOM predictor values")

# Transformation and scaling 
print("Transform the data  ... ")

# DEC predictors with abs(skew) > 1.
dec_skew <- apply(dec_dat, 2, skewness, na.rm = TRUE)
dec_skew[abs(dec_skew) > 1]

# Only one predictor from each season exceeds a skew of 2.0
# Those are fixed in place here. MakeMoreNormal() not required.
x <- dec_dat[, "dec_CS_bott_ave"]
y <- x^0.5
dec_dat[, "dec_CS_bott_ave"] <- y

x <- dec_dat[, "dec_CS_bott_max"]
y <- x^0.5
dec_dat[, "dec_CS_bott_max"] <- y

# Scaling December
x <- scale(dec_dat, center = T, scale = T)
dec_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist4 <- PlotHistos( dec_dat, "Transformed and scaled predictors" )
