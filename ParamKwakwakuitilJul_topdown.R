#---- Section ParamKwakwakuitilJul ----
# Pass 1: remove all cors > 0.9
showHighestPairs(jul_dat, 0.9)
#  jul_dat <- jul_dat[, !grepl("jul_SALT_surf_min", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_surf_min", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_bott_min", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_SALT_bott_max", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_surf_max", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_TEMP_bott_max", names(jul_dat))]

jul_dat <- jul_dat[, !grepl("jul_SALT_surf_ave", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_SALT_surf_min", names(jul_dat))]
jul_dat <- jul_dat[, !grepl("jul_SALT_bott_min", names(jul_dat))]

jul_dat <- jul_dat[, !grepl("jul_CS_bott_min", names(jul_dat))]

# Pass 2: remove all cors > 0.8
showHighestPairs(jul_dat, 0.8)

# Make histo of selected, unmodified predictors
hist1 <- PlotHistos( jul_dat, "Original FVCOM predictor values" )

# JUL predictors with abs(skew) > 1.
jul_skew <- apply(jul_dat, 2, skewness, na.rm = TRUE)
jul_skew[abs(jul_skew) > 2]

# Three predictors exceed a skew of 2.0
x <- jul_dat[, "jul_SALT_bott_ave"]
floor <- 20
y <- ifelse(x < floor, floor, x)
jul_dat[, "jul_SALT_bott_ave"] <- y

x <- jul_dat[, "jul_CS_bott_max"]
ceil <- 0.5
y <- ifelse(x > ceil, ceil, x)
jul_dat[, "jul_CS_bott_max"] <- y

x <- jul_dat[, "jul_CS_bott_ave"]
ceil <- 0.2
y <- ifelse(x > ceil, ceil, x)
jul_dat[, "jul_CS_bott_ave"] <- y

# Scaling July
x <- scale(jul_dat, center = T, scale = T)
jul_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist3 <- PlotHistos( jul_dat, "Transformed and scaled predictors" )





