#----  Parameter processing for Village Sea, Jul----

#Pass 1: remove all cors > 0.9
 showHighestPairs(jul_dat, 0.9)
jul_dat <- jul_dat[, !grepl("jul_SALT_bott_ave", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_bott_ave", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_SALT_surf_ave", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_surf_ave", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_surf_min", names(jul_dat))]

# Pass 2: remove all cors > 0.8
showHighestPairs(jul_dat, 0.8)
# Remove all July bottom currents
jul_dat <- jul_dat[, !grepl("jul_CS_bott", names(jul_dat))]

jul_dat <- jul_dat[, !grepl("jul_CS_surf_min", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_SALT_bott_min", names(jul_dat))]

# Pass 3: remove all cors > 0.5
# NOTE cor=0.7 and 0.6 fail cuz there is only 1 ...
showHighestPairs(jul_dat, 0.4)
# Cor btwn max bott temp and max bottom salt = -0.72 but both retain cuz chemistry.
names(jul_dat)

# Make histo of selected, unmodified predictors
hist1 <- PlotHistos(jul_dat, "Original FVCOM predictor values")

# JUL predictors with abs(skew) > 1.
jul_skew <- apply(jul_dat, 2, skewness, na.rm = TRUE)
jul_skew[abs(jul_skew) > 1]

# Only one predictor exceeds a skew of 2.0
x <- jul_dat[, "jul_SALT_surf_min"]
floor <- 21
y <- ifelse(x < floor, floor, x)
y <- y^3
jul_dat[, "jul_SALT_surf_min"] <- y

# Scaling July
x <- scale(jul_dat, center = T, scale = T)
dec_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist3 <- PlotHistos(jul_dat, "Transformed and scaled predictors")
