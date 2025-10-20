#---- Section ParamKwakwakuitilJul - Bottom up ----

# We start with the averages ... 
select <- jul_dat[, grepl("ave", names(jul_dat))]
names(select)

# Make histo of selected, unmodified predictors
hist1 <- PlotHistos( select, "Original FVCOM predictor values" )


# JUL predictors with abs(skew) > 1.
jul_skew <- apply(select, 2, skewness, na.rm = TRUE)
jul_skew[abs(jul_skew) > 1]

# TWO predictors exceed a skew of 1.0
x <- select[, "jul_SALT_bott_ave"]
floor <- 20
y <- ifelse(x < floor, floor, x)
y <- y^3
select[, "jul_SALT_bott_ave"] <- y


x <- select[, "jul_CS_bott_ave"]
ceil <- 0.2
y <- ifelse(x > ceil, ceil, x)
y <- y^0.5
select[, "jul_CS_bott_ave"] <- y

# Scaling July
x <- scale(select, center = T, scale = T)
jul_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist3 <- PlotHistos( jul_dat, "Transformed and scaled predictors" )


# Check cross-correlations
showHighestPairs(jul_dat, 0.5)

