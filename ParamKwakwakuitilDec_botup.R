#---- Section ParamKwakwakuitilDEC - Bottom up ----

# We start with the averages ... 
select <- dec_dat[, grepl("ave", names(dec_dat))]
names(select)

# Make histo of selected, unmodified predictors
hist2 <- PlotHistos( select, "Original FVCOM predictor values" )


# Predictors with abs(skew) > 1.
skew <- apply(select, 2, skewness, na.rm = TRUE)
skew[abs(skew) > 1]

# TWO predictors exceed a skew of 1.0
x <- select[, "dec_SALT_bott_ave"]
floor <- 20
y <- ifelse(x < floor, floor, x)
y <- y^3
select[, "dec_SALT_bott_ave"] <- y

x <- select[, "dec_CS_bott_ave"]
ceil <- 0.2
y <- ifelse(x > ceil, ceil, x)
y <- y^(1/2)
select[, "dec_CS_bott_ave"] <- y

# Scaling
x <- scale(select, center = T, scale = T)
dec_dat <- as.data.frame(x)

#Make histo of transformed and scaled predictors
hist4 <- PlotHistos( dec_dat, "Transformed and scaled predictors" )


# Check cross-correlations
showHighestPairs(dec_dat, 0.5)

