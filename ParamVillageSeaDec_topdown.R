# ---- Parameter processing Village Sea, Dec ----
  
# Pass 1: remove all cors > 0.9
showHighestPairs(dec_dat, 0.9)

# reduce chemistry correlations
dec_dat <- dec_dat[, !grepl("dec_SALT_bott_ave", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_TEMP_bott_ave", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_SALT_surf_ave", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_TEMP_surf_ave", names(dec_dat))]

# Pass 2: remove all cors > 0.8
showHighestPairs(dec_dat, 0.8)
# Remove all Dec bottom currents
dec_dat <- dec_dat[, !grepl("dec_CS_bott", names(dec_dat))]
dec_dat <- dec_dat[, !grepl("dec_CS_surf_min", names(dec_dat))]
# Remove min of temp correlation
dec_dat <- dec_dat[, !grepl("dec_TEMP_bott_min", names(dec_dat))]

# Pass 3: remove all cors > 0.7
showHighestPairs(dec_dat, 0.7)
# Remove min of salt ...
dec_dat <- dec_dat[, !grepl("dec_SALT_bott_min", names(dec_dat))]
# Remove max of temp...
dec_dat <- dec_dat[, !grepl("dec_TEMP_surf_max", names(dec_dat))]

# Pass 3: remove all cors > 0.6
showHighestPairs(dec_dat, 0.4)

#Make histo of selected, unmodified predictors
hist2 <- PlotHistos(dec_dat, "Original FVCOM predictor values")

# Transform for skew
x <- dec_dat[, "dec_SALT_surf_min"]
floor <- 21
y <- ifelse(x < floor, floor, x)
y <- y^8
dec_dat[, "dec_SALT_surf_min"] <- y

# Scaling December
x <- scale(dec_dat, center = T, scale = T)
dec_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist4 <- PlotHistos(dec_dat, "Transformed and scaled predictors")
