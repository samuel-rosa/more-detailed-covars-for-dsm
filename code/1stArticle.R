# DESCRIPTION ##################################################################
# Source code used to develop the first article of my PhD research project. The
# reference for the article is as follows:
# Samuel-Rosa, A.; Heuvelink, G. B. M.; Vasques, G. M. & Anjos, L. H. C. Do more
# detailed environmental covariates deliver more accurate soil maps?. Geoderma,
# v. 243-244, p. 214-227, 2015. doi:10.1016/j.geoderma.2014.12.017

# SETTINGS #####################################################################

# Clean workspace
rm(list = ls())
gc()

# Load packages
require(rgdal)
require(spgrass6)
require(raster)
require(gstat)
require(geoR)
require(rgeos)
require(pedometrics)
require(car)
require(caret)
require(MASS)
require(lattice)
require(latticeExtra)
require(grid)
require(gridExtra)
require(xtable)
require(Hmisc)
require(plotKML)
require(stringr)
require(plyr)
require(pbapply)
require(mail)

# Load auxiliary data
data(R_pal)
r_data <- "data/R/"
load(paste(r_data, "general.RData", sep = ""))

# Source user defined function
source(paste(r_code, "1stArticleHelper.R", sep = ""))

# Load and check calibration data (1:350)
cal_data <- paste(point_data, "labData.csv", sep = "")
cal_data <- read.table(cal_data, sep = ";", header = TRUE, dec = ".", 
                       na.strings = "na")
colnames(cal_data)
str(cal_data)
id <- c("sampleid", "longitude", "latitude", "CLAY", "ORCA", "ECEC")
id <- match(id, colnames(cal_data))
cal_data <- cal_data[1:350, id]
str(cal_data)
coordinates(cal_data) <- ~ longitude + latitude
proj4string(cal_data) <- sirgas2000
cal_data <- spTransform(cal_data, wgs1984utm22s)
plot(cal_data, pch = 20, cex = 0.5)

# Set path to results (figures and tables)
fig_dir <- "~/projects/dnos-sm-rs/res/fig/1stArticle/"
tab_dir <- "~/projects/dnos-sm-rs/res/tab/1stArticle/"
    
# Initiate GRASS GIS (64) section
initGRASS(gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS,
          location = "dnos-sm-rs", mapset = "predictions", 
          pid = Sys.getpid(), override = TRUE)
system("g.region rast=dnos.raster")
gmeta6()

# LOCATION OF THE STUDY AREA ###################################################

# Location of the study area in the Brazilian territory ------------------------
# This is the first figure of the article (FIG1a). We indicate some of the
# main Brazilian cities so that the reader can have a better picture of the 
# location of the study area compared to the size of the Brazilian territory.
brazil <- shapefile("~/dbGIS/dnos-sm-rs/brasil.shp")
bb <- bbox(brazil)
bb[1, 2] <- -34.7
brazil@bbox <- bb
pts <- data.frame(rbind(c(-53.790215, -29.657668), c(-47.887768, -15.788838),
                        c(-46.665337, -23.536445), c(-51.211134, -30.032484),
                        c(-59.992370, -3.080757)))
colnames(pts) <- c("long", "lat")
coordinates(pts) <- ~ long + lat
proj4string(pts) <- proj4string(brazil)
pts <- list("sp.points", pts, pch = 20, cex = 0.5, col = "black")
# set map colors
brazil$UF_05 <- as.factor(brazil$UF_05)
rs <- rep("lightgray", length(brazil$UF_05))
rs[which(brazil$UF_05 == "RS")] = "darkgray"
# prepare spplot
p <- spplot(brazil, zcol = "UF_05", aspect = "iso", col = "gray", 
            scales = list(draw = TRUE, 
                          x = list(at = seq(-70, -35, 5)),
                          y = list(at = seq(-30, 5, 5))),
            col.regions = colorRampPalette(rs)(27), colorkey = FALSE, 
            cex = 0.3, sp.layout = list(pts),
            par.settings = list(fontsize = list(text = 7, points = 5),
                                layout.widths = list(left.padding = 0, 
                                                     right.padding = 0), 
                                layout.heights = list(top.padding = 0,
                                                      bottom.padding = 0)),
            panel = function(x, y, ...) {
              panel.polygonsplot(x, y, ...)
              panel.abline(h = seq(-30, 0, 5), v = seq(-70, -40, 5),
                           col = "gray", lty = "dashed", lwd = 0.5)
              panel.text(x = -53.790215, y = -29.657668, "Santa Maria", pos = 2)
              panel.text(x = -47.887768, y = -15.788838, "Brasília", pos = 4)
              panel.text(x = -46.665337, y = -23.536445, "São Paulo", pos = 4)
              panel.text(x = -51.211134, y = -30.032484, "Porto Alegre", 
                         pos = 4)
              panel.text(x = -59.992370, y = -3.080757, "Manaus", pos = 4)
            }
)
p
# save image
dev.off()
pdf(file = paste(fig_dir, "FIG1a.pdf", sep = ""), width = 9/cm(1),
    height = 9/cm(1))
print(p)
dev.off()
rm(brazil, bb, pts, rs, p)
gc()

# Location of the calibration observations -------------------------------------
# This is the second figure of the article (FIG1b). We show the spatial 
# distribution of the calibration observations and of the drainage network. The
# last is included so that the reader can have a better picture of the 
# topography of the study area.
pol <- readVECT6("buffer_BASIN_10")
drain <- readVECT6("STREAM_10")
drain <- gIntersection(drain, pol, byid = TRUE)
p <- spplot(pol, zcol = "cat", col = "gray", fill = "lightgray",
            scales = list(draw = TRUE),
            colorkey = FALSE, aspect = "iso",
            xlim = c(bbox(cal_data)[1, 1] * 0.9998, 
                     bbox(cal_data)[1, 2] * 1.0005), 
            ylim = c(bbox(cal_data)[2, 1] * 0.99999,
                     bbox(cal_data)[2, 2] * 1.00001),
            par.settings = list(fontsize = list(text = 7, points = 5),
                                layout.widths = list(left.padding = 0, 
                                                     right.padding = 0), 
                                layout.heights = list(top.padding = 0,
                                                      bottom.padding = 0)),
            panel = function(x, y, ...) {
              panel.polygonsplot(x, y, ...)
              panel.points(x = coordinates(cal_data)[, 1],
                           y = coordinates(cal_data)[, 2],
                           pch = 20, cex = 0.5, col = "black")
              panel.abline(v = seq(227000, 232000, 1000), 
                           h = seq(6712000, 6722000, 1000),
                           col = "gray", lty = "dashed", lwd = 0.5)
            }) + layer(sp.lines(drain, col = "black", lty = "dashed", 
                                lwd = 0.3))
p
# save image
dev.off()
pdf(file = paste(fig_dir, "FIG1b.pdf", sep = ""), width = 9/cm(1),
    height = 9/cm(1))
print(p)
dev.off()
rm(pol, p, drain)
gc()

# convert pdf figures to png
pdf_file <- paste(fig_dir, "FIG1", letters[1:2], sep = "")
pdf2png(pdf_file)
rm(pdf_file)

# EXPLORATORY ANALYSIS #########################################################
# We check the soil data empirical distribution. Because it is severely skewed, 
# we transform it using the family of Box-Cox power transformations to achieve
# an empirical distribution closer to Gaussian. We also prepare histograms of
# frequency with the original and transformed variable.

# Box-Cox transformation
bc_lambda <- list(CLAY = NA, ORCA = NA, ECEC = NA)

# EXPLORATORY ANALYSIS - CLAY --------------------------------------------------
vari <- cal_data$CLAY
# Histogram with original variable
xlab <- expression(paste('CLAY (g ',kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", nint = 20, xlab = xlab, BoxCox = FALSE,
              col = c("lightgray", "black"), lty = "dashed",
              stats = FALSE, lwd = c(0.001, 0.5),
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(fig_dir, "FIG2a.pdf", sep = ""), width = 6.3/cm(1), 
    height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5), 
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Histogram with transformed variable
xlab  <-  expression(paste('Box-Cox CLAY (g ',kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)), lwd = c(0.001, 0.5))
dev.off()
pdf(file = paste(fig_dir, "FIG2d.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Transformation
lambda <- powerTransform(vari)$lambda
bc_lambda$CLAY <- lambda
cal_data$CLAY_BC <- bcPower(cal_data$CLAY, lambda)
xyplot(CLAY_BC ~ CLAY, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "Clay content")
rm(vari, lambda)
gc()

# EXPLORATORY ANALYSIS - ORCA --------------------------------------------------
vari <- cal_data$ORCA

# Histogram with original variable
xlab  <-  expression(paste('SOC (g ',kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", nint = 20, xlab = xlab, BoxCox = FALSE,
              col = c("lightgray", "black"), lty = "dashed",
              stats = FALSE, lwd = c(0.001, 0.5),
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(fig_dir, "FIG2b.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Histogram with transformed variable
xlab  <-  expression(paste('Box-Cox SOC (g ',kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)), lwd = c(0.001, 0.5))
dev.off()
pdf(file = paste(fig_dir, "FIG2e.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Transformation
lambda <- 0
bc_lambda$ORCA <- lambda
cal_data$ORCA_BC <- bcPower(cal_data$ORCA, lambda)
xyplot(ORCA_BC ~ ORCA, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "Carbon content")
rm(vari, lambda)
gc()

# EXPLORATORY ANALYSIS - ECEC --------------------------------------------------
vari <- cal_data$ECEC
# Histogram with original variable
xlab  <-  expression(paste('ECEC (',mmol, ' ', kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", nint = 20, xlab = xlab, BoxCox = FALSE,
              col = c("lightgray", "black"), lty = "dashed",
              stats = FALSE, lwd = c(0.001, 0.5),
              scales = list(cex = c(1, 1)))
dev.off()
pdf(file = paste(fig_dir, "FIG2c.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Histogram with transformed variable
xlab  <-  expression(paste('Box-Cox ECEC (',mmol, ' ', kg^-1,')', sep = ''))
tmp <- plotHD(vari, HD = "over", xlab = xlab, BoxCox = TRUE, stats = FALSE,
              scales = list(cex = c(1, 1)), lwd = c(0.001, 0.5))
dev.off()
pdf(file = paste(fig_dir, "FIG2f.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
print(tmp)
dev.off()
rm(tmp, xlab)
gc()
# Transformation
lambda <- 0
bc_lambda$ECEC <- lambda
cal_data$ECEC_BC <- bcPower(cal_data$ECEC, lambda)
xyplot(ECEC_BC ~ ECEC, data = cal_data@data, xlab = "original scale",
       ylab = "Box-Cox transformed", main = "ECEC")
rm(vari, lambda)
gc()

# Convert pdf figures to png
pdf_file <- paste(fig_dir, "FIG2", letters[1:6], sep = "")
pdf2png(pdf_file)
rm(pdf_file)

# COVARIATES ###################################################################
# We use two covariate databases that differ in their level of detail. Here we 
# define the covariates included in each database, prepare the sample points
# in the GRASS GIS database, and sample the raster files of each covariate 
# which are also stored in the GRASS GIS database. We do all the processng steps
# in GRASS GIS because loading all rasters in R would consume all the memory 
# available.

# COVARIATES - Database 1 ------------------------------------------------------
# These are the LESS detailed environmental covariates.
soil1 <- c("SOIL_100b", "SOIL_100c",
           "SOIL_100d", "SOIL_100e", "SOIL_100f")
soil1 <- paste(soil1, collapse = " + ")
land1 <- c("LU1980a", "LU1980b")
land1 <- paste(land1, collapse = " + ")
geo1 <- c("GEO_50a", "GEO_50b", "GEO_50c")
geo1 <- paste(geo1, collapse = " + ")
sat1 <- c("BLUE_30", "GREEN_30", "RED_30", "NIR_30a",
          "NIR_30b", "MIR_30", "NDVI_30", "SAVI_30")
sat1 <- paste(sat1, collapse = " + ")
dem1 <- c("ELEV_90", "SLP_90_3", "SLP_90_7", "SLP_90_15", "SLP_90_31", 
          "SLP_90_63", "SLP_90_127", "SLP_90_255", "TPI_90_3", "TPI_90_7",
          "TPI_90_15", "TPI_90_31", "TPI_90_63", "TPI_90_127", "TPI_90_255",
          "NOR_90_3", "NOR_90_7", "NOR_90_15", "NOR_90_31", "NOR_90_63",
          "NOR_90_127", "NOR_90_255", "TWI_90_3", "TWI_90_7", "TWI_90_15",
          "TWI_90_31", "TWI_90_63", "TWI_90_127", "TWI_90_255", "SPI_90_3",
          "SPI_90_7", "SPI_90_15", "SPI_90_31", "SPI_90_63", "SPI_90_127",
          "SPI_90_255")
dem1 <- paste(dem1, collapse = " + ")

# COVARIATES - Database 2 ------------------------------------------------------
# These are the MORE detailed environmental covariates.
soil2 <- c("SOIL_25a", "SOIL_25b", "SOIL_25c", "SOIL_25d", 
           "SOIL_25h", "SOIL_25i", "SOIL_25j")
soil2 <- paste(soil2, collapse = " + ")
land2 <- c("LU2009a", "LU2009b", "LU2009c", "LU2009d", "LUdiff")
land2 <- paste(land2, collapse = " + ")
geo2 <- c("GEO_25a", "GEO_25b", "GEO_25c", "GEO_25d")
geo2 <- paste(geo2, collapse = " + ")
sat2 <- c("BLUE_5", "GREEN_5", "RED_5", "EDGE_5", "NIR_5", "NDVI_5a",
          "NDVI_5b", "SAVI_5a", "SAVI_5b")
sat2 <- paste(sat2, collapse = " + ")
dem2 <- c("ELEV_10", "SLP_10_3", "SLP_10_7", "SLP_10_15", "SLP_10_31", 
          "SLP_10_63", "SLP_10_127", "SLP_10_255", "TPI_10_3", "TPI_10_7",
          "TPI_10_15", "TPI_10_31", "TPI_10_63", "TPI_10_127", "TPI_10_255",
          "NOR_10_3", "NOR_10_7", "NOR_10_15", "NOR_10_31", "NOR_10_63", 
          "NOR_10_127", "NOR_10_255", "TWI_10_3", "TWI_10_7", "TWI_10_15", 
          "TWI_10_31", "TWI_10_63", "TWI_10_127", "TWI_10_255", "SPI_10_3", 
          "SPI_10_7", "SPI_10_15", "SPI_10_31", "SPI_10_63", "SPI_10_127", 
          "SPI_10_255")
dem2 <- paste(dem2, collapse = " + ")

# COVARIATES - Sample rasters --------------------------------------------------
system("g.region rast=dnos.raster")
system("g.remove MASK")
system("r.mask dnos.raster")

# Import calibration points into GRASS
system("g.remove vect=calibration") # run twice
pts <- data.frame(coordinates(cal_data), cal_data$sampleid)
coordinates(pts) <- ~ longitude + latitude
proj4string(pts) <- wgs1984utm22s
colnames(pts@data) <- "sampleid"
spgrass6::writeVECT6(pts, "calibration", v.in.ogr_flags = "o")
rm(pts)

# Setup database of calibration points
system("v.info -c calibration")
cols_int <- paste(soil1, land1, geo1, soil2, land2, geo2, sep = " + ")
cols_int <- stringr::str_replace_all(cols_int, "[+]", "INT,")
cols_int <- paste(cols_int, " INT", sep = "")
cols_double <- paste(sat1, dem1,  sat2, dem2, sep = " + ")
cols_double <- stringr::str_replace_all(cols_double, "[+]", " DOUBLE PRECISION,")
cols_double <- paste(cols_double, " DOUBLE PRECISION", sep = "")
cols <- paste(cols_int, cols_double, sep = ", ")
cmd <- paste("v.db.addcol map=calibration columns='", cols, "'", sep = "")
system(cmd)
system("v.info -c calibration")
rm(cols_int, cols_double, cols)

# Sample rasters
maps <- paste(soil1, land1, geo1, soil2, land2, geo2, sat1, dem1, 
              sat2, dem2, sep = " + ")
maps <- stringr::str_replace_all(maps, "[ ]", "")
maps <- c(unlist(stringr::str_split(maps, "[+]")))
column <- maps
cmd <- paste("v.what.rast vect=calibration raster=", maps, " column=", 
             column, sep = "")
lapply(cmd, system)
rm(maps, column)

# Read calibration data into R
tmp <- spgrass6::readVECT6(vname = "calibration")
str(tmp)
tmp@data <- tmp@data[, - 1]
tmp$sampleid <- as.character(tmp$sampleid)
str(tmp)
which(c(cal_data$sampleid == tmp$sampleid) == FALSE)
sp::proj4string(tmp) <- wgs1984utm22s
tmp@data <- plyr::join(cal_data@data, tmp@data, by = "sampleid")
cal_data <- tmp
str(cal_data)
rm(tmp)

# CANDIDATE MODELS #############################################################

# CANDIDATE MODELS - Combinations ----------------------------------------------
# We start creating an object called 'combs' with all possible combinations of
# environmental covariates. The object is a list with two major items:
# - 'main' has the identification of the covariates;
# - 'num' is used for plotting purposes.
# The object also has two items used for a sensitivity analysis ('base' and 
# 'fine'). The sensitivity analysis shows the relative effect of excluding one 
# environmental covariate at a time from the candidate models.
combs <- list()
combs$main <- expand.grid(c("soil1", "soil2"), c("land1", "land2"),
                          c("geo1", "geo2"), c("sat1", "sat2"),
                          c("dem1", "dem2"), stringsAsFactors = FALSE)
combs$main <- split(combs$main, seq(1, nrow(combs$main), 1))
combs$num <- expand.grid(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2))
colnames(combs$num) <- c("soil", "land", "geo", "sat", "dem")
combs$base <- list()
for (i in 1:length(combs$main[[1]])) combs$base[[i]] <- combs$main[[1]][-i]
combs$fine <- list()
for (i in 1:length(combs$main[[32]])) combs$fine[[i]] <- combs$main[[32]][-i]
str(combs, 1)

# CANDIDATE MODELS - Predictors ------------------------------------------------
# Now we create a list with the predictor variables derived from the 
# environmental covariates.
preds <- list()
preds$main <- lapply(combs$main, function (X) parse(text = X))
preds$main <- lapply(preds$main, function (X) lapply(X, eval))
preds$main <- lapply(preds$main, function (X) paste(unlist(X), collapse = " + "))
preds$base <- lapply(combs$base, function (X) parse(text = X))
preds$base <- lapply(preds$base, function (X) lapply(X, eval))
preds$base <- lapply(preds$base, function (X) paste(unlist(X), collapse = " + "))
preds$fine <- lapply(combs$fine, function (X) parse(text = X))
preds$fine <- lapply(preds$fine, function (X) lapply(X, eval))
preds$fine <- lapply(preds$fine, function (X) paste(unlist(X), collapse = " + "))
str(preds, 1)

# CANDIDATE MODELS - Formulas --------------------------------------------------
# Now we use the lists created above to define a whole set of formulas that 
# define the candidate models for each soil property. This is done for the 
# sensitivity analysis too. We set 'y' as the response variable for all 
# formulas. Bellow, when we fit the models, we simply update the formulas using
# the respective soil property.
forms <- list()
forms$main <- lapply("y", function (X) paste(X, " ~ ", preds$main))
forms$main <- lapply(unlist(forms$main), as.formula)
forms$base <- lapply("y", function (X) paste(X, " ~ ", preds$base))
forms$base <- lapply(unlist(forms$base), as.formula)
forms$fine <- lapply("y", function (X) paste(X, " ~ ", preds$fine))
forms$fine <- lapply(unlist(forms$fine), as.formula)
str(forms, 1)

# LINEAR MODEL #################################################################
# Here we fit linear models using ordinary least squares to model the 
# deterministic component of the spatial variation. We use several strategies 
# to select the predictor variables so that we can evaluate the sensitivity of
# our results to the chosen method. The methods are:
# - use all predictors;
# - selection using the VIF;
# - selection using the VIF and stepwise AIC;
# - selection using the VIV and forward AIC;
# - selection using the VIV and backward AIC;
# The final linear models are selected using the VIF and stepwise AIC.

# LINEAR MODEL - CLAY ----------------------------------------------------------

# Update formulas
formula <- lapply(forms$main, update, CLAY_BC ~ .)
data <- cal_data@data

# Build linear model series
clay_full <- buildMS(formula, data)
clay_vif <- buildMS(formula, data, vif = TRUE)
clay_both <- buildMS(formula, data, vif = TRUE, aic = TRUE)
clay_for <- buildMS(formula, data, vif = TRUE, aic = TRUE,
                    aic.direction = "forward")
clay_back <- buildMS(formula, data, vif = TRUE, aic = TRUE, 
                     aic.direction = "backward")

# Get the statistics of the linear model series
clay_full_stats <- statsMS(clay_full, combs$num, "rmse")
clay_vif_stats <- statsMS(clay_vif, combs$num, "rmse")
clay_both_stats <- statsMS(clay_both, combs$num, "rmse")
clay_for_stats <- statsMS(clay_for, combs$num, "rmse")
clay_back_stats <- statsMS(clay_back, combs$num, "rmse")

# Save all linear model series plots
# This is to evaluate the behaviour of the results regarding the variable 
# selection method.
grid <- c(2:6)
line <- "ADJ_r2"
ind <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(clay_full_stats, grid, line, ind, color = color, 
                 main = "full model")
b_plot <- plotMS(clay_vif_stats, grid, line, ind, color = color, 
                 main = "VIF selection")
c_plot <- plotMS(clay_for_stats, grid, line, ind, color = color, 
                 main = "forward selection")
d_plot <- plotMS(clay_back_stats, grid, line, ind, color = color, 
                 main = "backward selection")
e_plot <- plotMS(clay_both_stats, grid, line, ind, color = color,
                 main = "stepwise selection")
dev.off()
pdf(file = paste(fig_dir, "clay_model_series_plot.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)

# Save the linear model series plot (pdf and png)
# We use the linear models calibrated using the stepwise variable selection
dev.off()
pdf(file = paste(fig_dir, "FIG5a.pdf", sep = ""),
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(clay_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"),
       xlab = "CLAY model ranking", scales = list(cex = c(1, 1)))
dev.off()
pdf2png(paste(fig_dir, "FIG5a", sep = ""))

# Check the effect of the number of observations
# This item has not been included in the published paper.
# We use three sample sizes (100, 200, and 300) and 50 iterations to calibrate
# linear models with the whole set of predictor variables. The 50 model
# series plots are used to produce an animated gif. The change observed in the
# animated gif among the model series plots is used as an indicator of the
# sensitivity of the results to the number of calibration observations.
# Change 'n' to set a different sample size.
n <- 300
for (i in 1:50) {
  data <- cal_data@data[sample(c(1:350), size = n), ]
  tmp <- buildMS(formula, data)
  tmp <- statsMS(tmp, combs$num, "rmse")
  color <- c("lightyellow", "palegreen")
  p <- plotMS(tmp, grid = c(2:6), line = "ADJ_r2", ind = 2, color = color)
  jpeg(paste("/tmp/plot", i, ".jpg", sep = ""))
  print(p)
  dev.off()
}
system("convert /tmp/*.jpg -delay 10 -loop 1 /tmp/movie.gif")
system("eog /tmp/movie.gif")

# Get base and best models
clay_sel <- list()
clay_sel$poor_lm <- clay_both[head(clay_both_stats, 1)$id][[1]]
clay_sel$base_lm <- clay_both[[1]]
clay_sel$fine_lm <- clay_both[[32]]
clay_sel$best_lm <- clay_both[tail(clay_both_stats, 1)$id][[1]]

#clay_base_lm <- clay_both[[1]]
#clay_best_lm <- rev(clay_both[tail(clay_both_stats, 1)$id])[[1]]

# Analysis of the residuals

# BASE MODEL
data <- cal_data
model <- clay_sel$base_lm
res <- residuals(model)
lambda <- bc_lambda$CLAY
par(mfrow = c(2, 3))
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])

# BEST MODEL
data  <- cal_data
model <- clay_sel$best_lm
res <- residuals(model)
lambda <- bc_lambda$clay
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])
dev.off()

# LINEAR MODEL - ORCA ----------------------------------------------------------

# Update formulas
formula <- lapply(forms$main, update, ORCA_BC ~ .)
data <- cal_data@data

# Build linear model series
orca_full <- buildMS(formula, data)
orca_vif <- buildMS(formula, data, vif = TRUE)
orca_both <- buildMS(formula, data, vif = TRUE, 
                     aic = TRUE, aic.direction = "both")
orca_for <- buildMS(formula, data, vif = TRUE, 
                    aic = TRUE, aic.direction = "forward")
orca_back <- buildMS(formula, data, vif = TRUE, 
                     aic = TRUE, aic.direction = "backward")

# Get the statistics of the linear model series
orca_full_stats <- statsMS(orca_full, combs$num, "rmse")
orca_vif_stats <- statsMS(orca_vif, combs$num, "rmse")
orca_both_stats <- statsMS(orca_both, combs$num, "rmse")
orca_for_stats <- statsMS(orca_for, combs$num, "rmse")
orca_back_stats <- statsMS(orca_back, combs$num, "rmse")

# Save all linear model series plots
# This is to evaluate the behaviour of the results regarding the variable 
# selection method.
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(orca_full_stats, grid, line, ind, color = color,
                 main = "full model")
b_plot <- plotMS(orca_vif_stats, grid, line, ind, color = color, 
                 main = "VIF selection")
c_plot <- plotMS(orca_for_stats, grid, line, ind, color = color, 
                 main = "forward selection")
d_plot <- plotMS(orca_back_stats, grid, line, ind, color = color, 
                 main = "backward selection")
e_plot <- plotMS(orca_both_stats, grid, line, ind, color = color, 
                 main = "stepwise selection")
dev.off()
pdf(file = paste(fig_dir, "orca_model_series_plot.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)
gc()

# Save the linear model series plot (pdf and png)
# We use the linear models calibrated using the stepwise variable selection
# A figure has been prepared for the presentation of Dr. Maria de Lourdes
# Mendonça Santos at the Brazilian Congress of Soil Science of 2015.
dev.off()
# pdf(file = "/home/alessandro/fig01.pdf", 
pdf(file = paste(fig_dir, "FIG5b.pdf", sep = ""),
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(orca_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2,
       # color = c("whitesmoke", "turquoise"),
       # ylab = c(expression(paste(R^2, " ajustado", sep = "")), "Covariáveis"),
       # xlab = paste("Posição dos modelos de predição do carbono orgânico do ",
                    # "solo\n\n", "Covariáveis: land - uso da terra, ",
                    # "geo - geologia, dem - atributos de terreno, soil - mapa ",
                    # "pedológico, sat - imagens de satélite\n",
                    # "Cores: cinza - modelos com covariáveis MENOS detalhadas,",
                    # " turquesa - modelos com covariáveis MAIS detalhadas\n",
                    # "Triângulos: posição média dos modelos com covariáveis ",
                    # "MAIS detalhadas", sep = ""),
       # scales = list(cex = c(1, 1)))
       color = c("lightyellow", "palegreen"),
       xlab = "SOC model ranking", scales = list(cex = c(1, 1)))
dev.off()
pdf2png(paste(fig_dir, "FIG5b", sep = ""))
# pdf2png("/home/alessandro/fig01")

# Check the effect of the number of observations
# This item has not been included in the published paper.
# We use three sample sizes (100, 200, and 300) and 50 iterations to calibrate
# linear models with the whole set of predictor variables. The 50 model
# series plots are used to produce an animated gif. The change observed in the
# animated gif among the model series plots is used as an indicator of the
# sensitivity of the results to the number of calibration observations.
# Change 'n' to set a different sample size.
n <- 300
for (i in 1:50) {
  data <- cal_data@data[sample(c(1:350), size = n), ]
  tmp <- buildMS(formula, data)
  tmp <- statsMS(tmp, combs$num, "rmse")
  color <- c("lightyellow", "palegreen")
  p <- plotMS(tmp, grid = c(2:6), line = "ADJ_r2", ind = 2, color = color)
  jpeg(paste("/tmp/plot", i, ".jpg", sep = ""))
  print(p)
  dev.off()
}
system("convert /tmp/*.jpg -delay 10 -loop 1 /tmp/movie.gif")
system("eog /tmp/movie.gif")

# Get base and best models
orca_sel <- list()
orca_sel$poor_lm <- orca_both[head(orca_both_stats, 1)$id][[1]]
orca_sel$base_lm <- orca_both[[1]]
orca_sel$fine_lm <- orca_both[[32]]
orca_sel$best_lm <- orca_both[tail(orca_both_stats, 1)$id][[1]]

# carbon_base_lm <- carbon_both[[1]]
# carbon_best_lm <- rev(carbon_both[tail(carbon_both_stats, 1)$id])[[1]]

# Analysis of the residuals

# BASE MODEL
data <- cal_data
model <- orca_sel$base_lm
res <- residuals(model)
lambda <- bc_lambda$ORCA
par(mfrow = c(2, 3))
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])

# BEST MODEL
data <- cal_data
model <- orca_sel$best_lm
res <- residuals(model)
lambda <- bc_lambda$carbon
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])
dev.off()

# LINEAR MODEL - ECEC ----------------------------------------------------------

# Update formulas
formula <- lapply(forms$main, update, ECEC_BC ~ .)
data <- cal_data@data

# Build linear model series
ecec_full <- buildMS(formula, data)
ecec_vif <- buildMS(formula, data, vif = TRUE)
ecec_both <- buildMS(formula, data, vif = TRUE, aic = TRUE, 
                     aic.direction = "both")
ecec_for <- buildMS(formula, data, vif = TRUE, aic = TRUE, 
                    aic.direction = "forward")
ecec_back <- buildMS(formula, data, vif = TRUE, aic = TRUE, 
                     aic.direction = "backward")

# Get the statistics of the linear model series
ecec_full_stats <- statsMS(ecec_full, combs$num, "rmse")
ecec_vif_stats <- statsMS(ecec_vif, combs$num, "rmse")
ecec_both_stats <- statsMS(ecec_both, combs$num, "rmse")
ecec_for_stats <- statsMS(ecec_for, combs$num, "rmse")
ecec_back_stats <- statsMS(ecec_back, combs$num, "rmse")

# Save all linear model series plots
# This is to evaluate the behaviour of the results regarding the variable 
# selection method.
grid <- c(2:6)
line <- "ADJ_r2"
ind  <- 2
color <- c("lightyellow", "palegreen")
a_plot <- plotMS(ecec_full_stats, grid, line, ind, color = color, 
                 main = "full model")
b_plot <- plotMS(ecec_vif_stats, grid, line, ind, color = color, 
                 main = "VIF selection")
c_plot <- plotMS(ecec_for_stats, grid, line, ind, color = color, 
                 main = "forward selection")
d_plot <- plotMS(ecec_back_stats, grid, line, ind, color = color, 
                 main = "backward selection")
e_plot <- plotMS(ecec_both_stats, grid, line, ind, color = color, 
                 main = "stepwise selection")
dev.off()
pdf(file = paste(fig_dir, "ecec_model_series_plot.pdf", sep = ""),
    width = 7, height = 15)
trellis.par.set(fontsize = list(text = 8, points = 6))
grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot, ncol = 1)
dev.off()
rm(grid, line, ind, color, a_plot, b_plot, c_plot, d_plot, e_plot)
gc()

# Save the linear model series plot (pdf and png)
# We use the linear models calibrated using the stepwise variable selection
dev.off()
pdf(file = paste(fig_dir, "FIG5c.pdf", sep = ""), 
    width = 19/cm(1), height = 8/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plotMS(ecec_both_stats, grid = c(2:6), line = "ADJ_r2", ind = 2, 
       color = c("lightyellow", "palegreen"),
       xlab = "ECEC model ranking", scales = list(cex = c(1, 1)))
dev.off()
pdf2png(paste(fig_dir, "FIG5c", sep = ""))

# Check the effect of the number of observations
# This item has not been included in the published paper.
# We use three sample sizes (100, 200, and 300) and 50 iterations to calibrate
# linear models with the whole set of predictor variables. The 50 model
# series plots are used to produce an animated gif. The change observed in the
# animated gif among the model series plots is used as an indicator of the
# sensitivity of the results to the number of calibration observations.
# Change 'n' to set a different sample size.
n <- 100
for (i in 1:50) {
  data <- cal_data@data[sample(c(1:350), size = n), ]
  tmp <- buildMS(formula, data)
  tmp <- statsMS(tmp, combs$num, "rmse")
  color <- c("lightyellow", "palegreen")
  p <- plotMS(tmp, grid = c(2:6), line = "ADJ_r2", ind = 2, color = color)
  jpeg(paste("/tmp/plot", i, ".jpg", sep = ""))
  print(p)
  dev.off()
}
system("convert /tmp/*.jpg -delay 10 -loop 1 /tmp/movie.gif")
system("eog /tmp/movie.gif")

# Get base and best models
ecec_sel <- list()
ecec_sel$poor_lm <- ecec_both[head(ecec_both_stats, 1)$id][[1]]
ecec_sel$base_lm <- ecec_both[[1]]
ecec_sel$fine_lm <- ecec_both[[32]]
ecec_sel$best_lm <- ecec_both[tail(ecec_both_stats, 1)$id][[1]]

# ecec_base_lm <- ecec_both[[1]]
# ecec_best_lm <- rev(ecec_both[tail(ecec_both_stats, 1)$id])[[1]]

# Analysis of the residuals

# BASE MODEL
data <- cal_data
model <- ecec_sel$base_lm
res <- residuals(model)
par(mfrow = c(2, 3))
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])

# BEST MODEL
data <- cal_data
model <- ecec_sel$best_lm
res <- residuals(model)
plot(model, which = c(1:6))
plotESDA(res, lon = coordinates(data)[, 1], lat = coordinates(data)[, 2])
dev.off()
rm(model, res)

# SENSITIVITY ANALYSIS #########################################################
# We evaluate here the change (increase or decrease) of the importance of each
# environmental covariate on explaining the variance when the more detailed
# version was used.

# Effect of dropping one environmental covariate -------------------------------
data <- cal_data@data
drop <- list()

# clay
formula <- lapply(forms$base, update, CLAY_BC ~ .)
drop$clay$base_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
formula <- lapply(forms$fine, update, CLAY_BC ~ .)
drop$clay$fine_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
drop$clay$base_r2 <- statsMS(drop$clay$base_lm)
drop$clay$fine_r2 <- statsMS(drop$clay$fine_lm)
drop$clay$base_dr2 <- deltaR2(clay_both_stats[clay_both_stats$id == 1, ],
                              drop$clay$base_r2)
drop$clay$fine_dr2 <- deltaR2(clay_both_stats[clay_both_stats$id == 32, ],
                              drop$clay$fine_r2)

# orca
formula <- lapply(forms$base, update, ORCA_BC ~ .)
drop$orca$base_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
formula <- lapply(forms$fine, update, ORCA_BC ~ .)
drop$orca$fine_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
drop$orca$base_r2 <- statsMS(drop$orca$base_lm)
drop$orca$fine_r2 <- statsMS(drop$orca$fine_lm)
drop$orca$base_dr2 <- deltaR2(orca_both_stats[orca_both_stats$id == 1, ],
                              drop$orca$base_r2)
drop$orca$fine_dr2 <- deltaR2(orca_both_stats[orca_both_stats$id == 32, ],
                              drop$orca$fine_r2)

# ecec
formula <- lapply(forms$base, update, ECEC_BC ~ .)
drop$ecec$base_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
formula <- lapply(forms$fine, update, ECEC_BC ~ .)
drop$ecec$fine_lm <- buildMS(formula, data, vif = TRUE, aic = TRUE)
rm(formula)
drop$ecec$base_r2 <- statsMS(drop$ecec$base_lm)
drop$ecec$fine_r2 <- statsMS(drop$ecec$fine_lm)
drop$ecec$base_dr2 <- deltaR2(ecec_both_stats[ecec_both_stats$id == 1, ],
                              drop$ecec$base_r2)
drop$ecec$fine_dr2 <- deltaR2(ecec_both_stats[ecec_both_stats$id == 32, ],
                              drop$ecec$fine_r2)
rm(data)

# Save LaTeX table
Covariate <- data.frame(drop$clay$base_dr2$ADJ_r2, drop$clay$fine_dr2$ADJ_r2,
                        drop$orca$base_dr2$ADJ_r2, drop$orca$fine_dr2$ADJ_r2,
                        drop$ecec$base_dr2$ADJ_r2, drop$ecec$fine_dr2$ADJ_r2)
Covariate <- round(Covariate, 3)
rownames(Covariate) <- c("\\texttt{soil}", "\\texttt{land}",
                      "\\texttt{geo}", "\\texttt{sat}", "\\texttt{dem}")
colnames(Covariate) <- rep(c("Less", "More"), 3)
long_cap <- paste("The importance of each environmental covariate$^a$ ",
                  "($\\Delta{R}^{2}_{adj}{}^b$) in the models calibrated ",
                  "with their less and more spatially detailed version.",
                  sep = "")
foot <- paste("${}^a$ Covariate: \\texttt{soil} - soil map, \\texttt{land} - ",
              "land use map, \\texttt{geo} - geologic map, \\texttt{sat} - ",
              "satellite image, and \\texttt{dem} - digital elevation model.",
              " ${}^b$ $\\Delta{R}^{2}_{adj} = {R}^{2}_{adj}{}_{q=5} - ",
              "R^{2}_{adj}{}_{q=5-1}$, where $q$ is the number of covariates ",
              "included in the model. Negative values result from adjusting ", 
              "the $R^{2}$ using the number of predictor variables initially ",
              "offered to enter the model instead of the reduced number of ",
              "predictor variables that entered the model.", sep = "")
file <- paste(tab_dir, "TAB4.tex", sep = "")
latex(Covariate, file = file, label = "tab:drop", table.env = TRUE, 
      longtable = FALSE, cgroup = c("CLAY","SOC", "ECEC"),
      n.cgroup = c(2, 2, 2), na.blank = TRUE, ctable = TRUE, caption = long_cap,
      where = "!h", insert.bottom = foot, rowlabel.just = "l", 
      cgroup.just = "l", cgroupTexCmd = NULL, rgroupTexCmd = NULL)
tmp <- readLines(file)
tmp <- gsub("pos=!h,]", "pos=!h,doinside={\\scriptsize\\setstretch{1.1}}]", tmp)
writeLines(tmp, file)
rm(Covariate, long_cap, foot, file)

# LINEAR MIXED MODEL ###########################################################
# We fit linear mixed models using the poor, base, fine, and best linear models
# fitter above. The variogram model fit to the empirical variogram is saved to
# compose the 6th figure of the article. Only the variograms of the base and 
# best linear mixed models are used.
clay_vario <- list()
orca_vario <- list()
ecec_vario <- list()

# LINEAR MIXED MODEL - clay ----------------------------------------------------
y <- "CLAY"
lambda <- bc_lambda$CLAY

## Poor linear mixed model
# Calculate and check the empirical variogram
model <- clay_sel$poor_lm
breaks <- seq(0, 6500, 100)
clay_vario$poor <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(clay_vario$poor, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.9, 1.0, 1.1), c(400, 500, 600)))
nugget <- c(0.3, 0.4, 0.5)
clay_sel$poor_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(clay_sel$poor_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Base linear mixed model

# Calculate and check the empirical variogram
model <- clay_sel$base_lm
breaks <- seq(0, 6500, 100)
clay_vario$base <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(clay_vario$base, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.9, 1.0, 1.1), c(400, 500, 600)))
nugget <- c(0.3, 0.4, 0.5)
clay_sel$base_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(clay_sel$base_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Fine linear mixed model

# Calculate and check the empirical variogram
model <- clay_sel$fine_lm
breaks <- seq(0, 6500, 100)
clay_vario$fine <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(clay_vario$fine, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(1.0, 1.1, 1.2), c(300, 400, 500)))
nugget <- c(0.1, 0.2, 0.3)
clay_sel$fine_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(clay_sel$fine_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Best linear mixed model

# Calculate and check the empirical variogram
model <- clay_sel$best_lm
breaks <- seq(0, 6500, 100)
clay_vario$best <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(clay_vario$best, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(1.0, 1.1, 1.2), c(300, 400, 500)))
nugget <- c(0.1, 0.2, 0.3)
clay_sel$best_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(clay_sel$best_lmm)
rm(model, breaks, nugget, tmp, ini.cov.pars)
gc()

# CLAY - plot experimental variograms and fitted linear mixed models -----------
xlim <- c(0, 3000)
ylim <- max(#clay_vario$poor$v[clay_vario$poor$u <= max(xlim)],
            clay_vario$base$v[clay_vario$base$u <= max(xlim)],
            #clay_vario$fine$v[clay_vario$fine$u <= max(xlim)],
            clay_vario$best$v[clay_vario$best$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 1))
#l1 <- linesREML(clay_sel$poor_lmm, add = FALSE)
l2 <- linesREML(clay_sel$base_lmm, add = FALSE)
#l3 <- linesREML(clay_sel$fine_lmm, add = FALSE)
l4 <- linesREML(clay_sel$best_lmm, add = FALSE)
# v1 <- xyplot(clay_poor_vario$v ~ clay_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, lty = 2, col = "black")
#              })
clay_v2 <- xyplot(clay_vario$base$v ~ clay_vario$base$u, ylim = ylim, pch = 20,
             scales = list(tick.number = 8), xlim = xlim, col =  "black",
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l2$x, y = l2$y, lty = 2, col = "black")
             })
# v3 <- xyplot(clay_fine_vario$v ~ clay_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, lty = 2, col = "black")
#              })
clay_v4 <- xyplot(clay_vario$best$v ~ clay_vario$best$u, ylim = ylim, pch = 20,
             scales = list(tick.number = 8), xlim = xlim, col =  "black", 
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               panel.lines(x = l4$x, y = l4$y, lty = 2, col = "black")
             })
ylab <- expression(paste("Semivariance [(g kg"^"-1", ")"^"2", "]", sep = ""))
# v1 <- update(c(v1, v2, v3, v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
#             scales = list(cex = c(0.9, 0.9)), layout = c(4, 1))
clay_v <- update(c(clay_v2, clay_v4), ylab = ylab, xlab = "Distance [m]", 
                 asp = 1, scales = list(cex = c(1, 1)), layout = c(2, 1))
# save plot
dev.off()
pdf(file = paste(fig_dir, "FIG6a.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(clay_v)
dev.off()
#rm(xlim, ylim, l1, l2, l3, l4, clay_v1, clay_v2, clay_v3, clay_v4, clay_v)
rm(xlim, ylim, l2, l4, clay_v2, clay_v4, clay_v)
gc()

# LINEAR MIXED MODEL - orca ----------------------------------------------------
y <- "ORCA"
lambda <- bc_lambda$ORCA

## Poor linear mixed model
# Calculate and check the empirical variogram
model <- orca_sel$poor_lm
breaks <- seq(0, 6500, 100)
orca_vario$poor <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(orca_vario$poor, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget <- c(0.10, 0.12, 0.14)
orca_sel$poor_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(orca_sel$poor_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Base linear mixed model

# Calculate and check the empirical variogram
model <- orca_sel$base_lm
breaks <- seq(0, 6500, 100)
orca_vario$base <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(orca_vario$base, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget <- c(0.10, 0.12, 0.14)
orca_sel$base_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(orca_sel$base_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Fine linear mixed model

# Calculate and check the empirical variogram
model <- orca_sel$fine_lm
breaks <- seq(0, 6500, 100)
orca_vario$fine <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(orca_vario$fine, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.12, 0.14, 0.16), c(500, 600, 700)))
nugget <- c(0.10, 0.12, 0.14)
orca_sel$fine_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(orca_sel$fine_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Best linear mixed model

# Calculate and check the empirical variogram
model <- orca_sel$best_lm
breaks <- seq(0, 6500, 100)
orca_vario$best <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(orca_vario$best, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.10, 0.12, 0.14), c(400, 500, 600)))
nugget <- c(0.10, 0.12, 0.14)
orca_sel$best_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(orca_sel$best_lmm)
rm(model, breaks, nugget, tmp, ini.cov.pars)
gc()

# ORCA - plot experimental variograms and fitted linear mixed models -----------
xlim <- c(0, 3000)
ylim <- max(#orca_vario$poor$v[orca_vario$poor$u <= max(xlim)],
            orca_vario$base$v[orca_vario$base$u <= max(xlim)],
            #orca_vario$fine$v[orca_vario$fine$u <= max(xlim)],
            orca_vario$best$v[orca_vario$best$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 1))
#l1 <- linesREML(orca_sel$poor_lmm, add = FALSE)
l2 <- linesREML(orca_sel$base_lmm, add = FALSE)
#l3 <- linesREML(orca_sel$fine_lmm, add = FALSE)
l4 <- linesREML(orca_sel$best_lmm, add = FALSE)
# v1 <- xyplot(orca_poor_vario$v ~ orca_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, lty = 2, col = "black")
#              })
orca_v2 <- xyplot(orca_vario$base$v ~ orca_vario$base$u, ylim = ylim, pch = 20,
                  scales = list(tick.number = 8), xlim = xlim, col =  "black",
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    panel.lines(x = l2$x, y = l2$y, lty = 2, col = "black")
                  })
# v3 <- xyplot(orca_fine_vario$v ~ orca_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, lty = 2, col = "black")
#              })
orca_v4 <- xyplot(orca_vario$best$v ~ orca_vario$best$u, ylim = ylim, pch = 20,
                  scales = list(tick.number = 8), xlim = xlim, col =  "black", 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    panel.lines(x = l4$x, y = l4$y, lty = 2, col = "black")
                  })
ylab <- expression(paste("Semivariance [(g kg"^"-1", ")"^"2", "]", sep = ""))
# v1 <- update(c(v1, v2, v3, v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
#             scales = list(cex = c(0.9, 0.9)), layout = c(4, 1))
orca_v <- update(c(orca_v2, orca_v4), ylab = ylab, xlab = "Distance [m]", 
                 asp = 1, scales = list(cex = c(1, 1)), layout = c(2, 1))
# save plot
dev.off()
pdf(file = paste(fig_dir, "FIG6b.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(orca_v)
dev.off()
#rm(xlim, ylim, l1, l2, l3, l4, orca_v1, orca_v2, orca_v3, orca_v4, orca_v)
rm(xlim, ylim, l2, l4, orca_v2, orca_v4, orca_v)
gc()

# LINEAR MIXED MODEL - ecec ----------------------------------------------------
y <- "ECEC"
lambda <- bc_lambda$ECEC

## Poor linear mixed model
# Calculate and check the empirical variogram
model <- ecec_sel$poor_lm
breaks <- seq(0, 6500, 100)
ecec_vario$poor <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(ecec_vario$poor, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.30, 0.32, 0.34), c(500, 600, 700)))
nugget <- c(0.20, 0.22, 0.24)
ecec_sel$poor_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(ecec_sel$poor_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Base linear mixed model

# Calculate and check the empirical variogram
model <- ecec_sel$base_lm
breaks <- seq(0, 6500, 100)
ecec_vario$base <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(ecec_vario$base, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.30, 0.32, 0.34), c(500, 600, 700)))
nugget <- c(0.20, 0.22, 0.24)
ecec_sel$base_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(ecec_sel$base_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Fine linear mixed model

# Calculate and check the empirical variogram
model <- ecec_sel$fine_lm
breaks <- seq(0, 6500, 100)
ecec_vario$fine <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(ecec_vario$fine, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.20, 0.24, 0.28), c(500, 600, 700)))
nugget <- c(0.16, 0.20, 0.24)
ecec_sel$fine_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(ecec_sel$fine_lmm)
rm(model, breaks, nugget, ini.cov.pars)
gc()

## Best linear mixed model

# Calculate and check the empirical variogram
model <- ecec_sel$best_lm
breaks <- seq(0, 6500, 100)
ecec_vario$best <- fitVariog(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             lambda = lambda, breaks = breaks)
plot(ecec_vario$best, col = "blue", pch = 20, type = "b")

# Estimate model parameters using REML
ini.cov.pars <- as.matrix(expand.grid(c(0.20, 0.24, 0.28), c(500, 600, 700)))
nugget <- c(0.16, 0.20, 0.24)
ecec_sel$best_lmm <- fitREML(y = y, model = model, 
                             data = as.data.frame(cal_data), 
                             ini.cov.pars = ini.cov.pars,
                             nugget = nugget, lambda = lambda)
summary(ecec_sel$best_lmm)
rm(model, breaks, nugget, tmp, ini.cov.pars)
gc()

# ECEC - plot experimental variograms and fitted linear mixed models -----------
xlim <- c(0, 3000)
ylim <- max(#ecec_vario$poor$v[ecec_vario$poor$u <= max(xlim)],
  ecec_vario$base$v[ecec_vario$base$u <= max(xlim)],
  #ecec_vario$fine$v[ecec_vario$fine$u <= max(xlim)],
  ecec_vario$best$v[ecec_vario$best$u <= max(xlim)]) * 1.1
ylim <- c(0, round(ylim, 1))
#l1 <- linesREML(ecec_sel$poor_lmm, add = FALSE)
l2 <- linesREML(ecec_sel$base_lmm, add = FALSE)
#l3 <- linesREML(ecec_sel$fine_lmm, add = FALSE)
l4 <- linesREML(ecec_sel$best_lmm, add = FALSE)
# v1 <- xyplot(ecec_poor_vario$v ~ ecec_poor_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l1$x, y = l1$y, lty = 2, col = "black")
#              })
ecec_v2 <- xyplot(ecec_vario$base$v ~ ecec_vario$base$u, ylim = ylim, pch = 20,
                  scales = list(tick.number = 8), xlim = xlim, col =  "black",
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    panel.lines(x = l2$x, y = l2$y, lty = 2, col = "black")
                  })
# v3 <- xyplot(ecec_fine_vario$v ~ ecec_fine_vario$u, ylim = ylim, pch = 20,
#              col =  "black", scales = list(tick.number = 8), xlim = xlim,
#              panel = function(x, y, ...) {
#                panel.xyplot(x, y, ...)
#                panel.lines(x = l3$x, y = l3$y, lty = 2, col = "black")
#              })
ecec_v4 <- xyplot(ecec_vario$best$v ~ ecec_vario$best$u, ylim = ylim, pch = 20,
                  scales = list(tick.number = 8), xlim = xlim, col =  "black", 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, ...)
                    panel.lines(x = l4$x, y = l4$y, lty = 2, col = "black")
                  })
ylab <- expression(paste("Semivariance [(mmol kg"^"-1", ")"^"2", "]", sep = ""))
# v1 <- update(c(v1, v2, v3, v4), ylab = ylab, xlab = "Distance [m]", asp = 1,
#             scales = list(cex = c(0.9, 0.9)), layout = c(4, 1))
ecec_v <- update(c(ecec_v2, ecec_v4), ylab = ylab, xlab = "Distance [m]", 
                 asp = 1, scales = list(cex = c(1, 1)), layout = c(2, 1))
# save plot
dev.off()
pdf(file = paste(fig_dir, "FIG6c.pdf", sep = ""), 
    width = 9/cm(1), height = 6/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(ecec_v)
dev.off()
#rm(xlim, ylim, l1, l2, l3, l4, ecec_v1, ecec_v2, ecec_v3, ecec_v4, ecec_v)
rm(xlim, ylim, l2, l4, ecec_v2, ecec_v4, ecec_v)
gc()

# Convert pdf figures to png
pdf_file <- paste(fig_dir, "FIG6", letters[1:3], sep = "")
pdf2png(pdf_file)
rm(pdf_file)

# SAVE - LINEAR MODELS #########################################################
# We save all data produced till now and remove large objects containing linear
# models calibrated using different covariate selection procedures.
save(list = ls(), file = paste(r_data, "1stArticlePartI.rda", sep = ""))
rm(clay_back, clay_back_stats, clay_both, clay_both_stats, clay_for, 
   clay_for_stats, clay_full, clay_full_stats, clay_vif, clay_vif_stats, drop,
   ecec_back, ecec_back_stats, ecec_both, ecec_both_stats, ecec_for, 
   ecec_for_stats, ecec_full, ecec_full_stats, ecec_vif, ecec_vif_stats,
   orca_back, orca_back_stats, orca_both, orca_both_stats, orca_for, 
   orca_for_stats, orca_full, orca_full_stats, orca_vif, orca_vif_stats)
gc()

# Load data to continue working
load(file = "data/R/1stArticlePartI.rda")
ls()
rm(clay_back, clay_back_stats, clay_both, clay_both_stats, clay_for, 
   clay_for_stats, clay_full, clay_full_stats, clay_vif, clay_vif_stats, drop,
   ecec_back, ecec_back_stats, ecec_both, ecec_both_stats, ecec_for, 
   ecec_for_stats, ecec_full, ecec_full_stats, ecec_vif, ecec_vif_stats,
   orca_back, orca_back_stats, orca_both, orca_both_stats, orca_for, 
   orca_for_stats, orca_full, orca_full_stats, orca_vif, orca_vif_stats)
gc()

# Save best linear mixed models for the second article
clay_lmm <- clay_sel$best_lmm
orca_lmm <- orca_sel$best_lmm
ecec_lmm <- ecec_sel$best_lmm
save(clay_lmm, orca_lmm, ecec_lmm, file = paste(r_data, "2stArticlePartI.rda"))
rm(clay_lmm, orca_lmm, ecec_lmm)
gc()

# LEAVE-ONE-OUT CROSS-VALIDATION ###############################################
# We check our models using leave-one-out cross-validation. The predictions are
# back-transformed using numerical simulation. In both cases (lm and lmm) all
# model parameters are reestimated at each cross-validation step, significantly 
# increasing the computation time. Because we used an authomated covariate
# selection, it would be more elegant to submit the whole set of covariates 
# to the authomated selection at each cross-validation step. However, we have 
# chosen not to do so because we saw above with the model series plots 
# that the results are quite stable independent of the covariate selection 
# method used. We stress that during the back-transformation
# of values predicted using the linear model we assume that the prediction error
# variance is estimated without bias.

# Check the needed number of realizations --------------------------------------
# We start checking how many realizations are needed to stabilize the variance
# when back-transforming the predictions from the Box-CoX space to the original
# soil data space. We use the same number of realizations for all cases.
#  The results show that we should use at least 20,000 realizations.

# CLAY
model <- clay_sel$base_lm
lambda <- bc_lambda$CLAY
clay_sel$base_lm_cv <- looCV(model = model, simul.back = FALSE)
invBoxCox(mean = clay_sel$base_lm_cv$pred, variance = clay_sel$base_lm_cv$pev,
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda)

# SOC
model <- orca_sel$base_lm
lambda <- bc_lambda$ORCA
orca_sel$base_lm_cv <- looCV(model = model, simul.back = FALSE)
invBoxCox(mean = orca_sel$base_lm_cv$pred, variance = orca_sel$base_lm_cv$pev,
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda)

# ECEC
model <- ecec_sel$base_lm
lambda <- bc_lambda$ECEC
ecec_sel$base_lm_cv <- looCV(model = model, simul.back = FALSE)
invBoxCox(mean = ecec_sel$base_lm_cv$pred, variance = ecec_sel$base_lm_cv$pev,
          lambda = lambda, profile = TRUE, simul.back = TRUE)
rm(model, lambda)

# CLAY - base linear model -----------------------------------------------------
model <- clay_sel$base_lm
lambda <- bc_lambda$CLAY
clay_sel$base_lm_cv <- looCV(model = model, back = TRUE, lambda = lambda,
                             original = cal_data$CLAY, simul.back = TRUE, 
                             n.sim = 20000)
rm(model, lambda)

# CLAY - best linear model -----------------------------------------------------
model <- clay_sel$best_lm
lambda <- bc_lambda$CLAY
clay_sel$best_lm_cv <- looCV(model = model, back = TRUE, lambda = lambda,
                             original = cal_data$CLAY, simul.back = TRUE, 
                             n.sim = 20000)
rm(model, lambda)

# CLAY - base linear mixed model -----------------------------------------------
model <- clay_sel$base_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- clay_sel$base_lmm
clay_sel$base_lmm_cv <- krigeCV(geodata = geodata, model = model, 
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)

# CLAY - best linear mixed model -----------------------------------------------
model <- clay_sel$best_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- clay_sel$best_lmm
clay_sel$best_lmm_cv <- krigeCV(geodata = geodata, model = model, 
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)

# CARBON - base linear model ---------------------------------------------------
model <- orca_sel$base_lm
lambda <- bc_lambda$ORCA
orca_sel$base_lm_cv <- looCV(model = model, back = TRUE, simul.back = TRUE, 
                             original = cal_data$ORCA, lambda = lambda, 
                             n.sim = 20000)
rm(model, lambda)

# CARBON - best linear model ---------------------------------------------------
model <- orca_sel$best_lm
lambda <- bc_lambda$ORCA
orca_sel$best_lm_cv <- looCV(model = model, back = TRUE, simul.back = TRUE, 
                             original = cal_data$ORCA, lambda = lambda,
                             n.sim = 20000)
rm(model, lambda)

# CARBON - base linear mixed model ---------------------------------------------
model <- orca_sel$base_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model<- orca_sel$base_lmm
orca_sel$base_lmm_cv <- krigeCV(geodata = geodata, model = model, 
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)

# CARBON - best linear mixed model ---------------------------------------------
model <- orca_sel$best_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model<- orca_sel$best_lmm
orca_sel$best_lmm_cv <- krigeCV(geodata = geodata, model = model, 
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)

# ECEC - base linear model -----------------------------------------------------
model <- ecec_sel$base_lm
lambda <- bc_lambda$ECEC
ecec_sel$base_lm_cv <- looCV(model = model, back = TRUE, simul.back = TRUE, 
                             original = cal_data$ECEC, lambda = lambda, 
                             n.sim = 20000)
rm(model, lambda)

# ECEC - best linear model -----------------------------------------------------
model <- ecec_sel$best_lm
lambda <- bc_lambda$ECEC
ecec_sel$best_lm_cv <- looCV(model = model, back = TRUE, simul.back = TRUE, 
                             lambda = lambda, original = cal_data$ECEC,
                             n.sim = 20000)
rm(model, lambda)

# ECEC - base linear mixed model -----------------------------------------------
model <- ecec_sel$base_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- ecec_sel$base_lmm
ecec_sel$base_lmm_cv <- krigeCV(model = model, geodata = geodata, back = TRUE,
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)
gc()

# ECEC - best linear mixed model -----------------------------------------------
model <- ecec_sel$best_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- ecec_sel$best_lmm
ecec_sel$best_lmm_cv <- krigeCV(model = model, geodata = geodata, back = TRUE,
                                reestimate = TRUE, output.reestimate = TRUE, 
                                n.sim = 20000)
rm(geodata, model)
gc()

# Tables with results ----------------------------------------------------------

# Cross-validation statistics (ME, MAE, RMSE, SRMSE, AVE)
obs <- data.frame(cal_data$CLAY, cal_data$CLAY, cal_data$CLAY, cal_data$CLAY, 
                  cal_data$ORCA, cal_data$ORCA, cal_data$ORCA, cal_data$ORCA,
                  cal_data$ECEC, cal_data$ECEC, cal_data$ECEC, cal_data$ECEC)
str(obs)               
pred <- data.frame(clay_sel$base_lm_cv$back.transformed$pred,
                   clay_sel$base_lmm_cv$back.transformed$pred,
                   clay_sel$best_lm_cv$back.transformed$pred,
                   clay_sel$best_lmm_cv$back.transformed$pred,
                   orca_sel$base_lm_cv$back.transformed$pred,
                   orca_sel$base_lmm_cv$back.transformed$pred,
                   orca_sel$best_lm_cv$back.transformed$pred,
                   orca_sel$best_lmm_cv$back.transformed$pred,
                   ecec_sel$base_lm_cv$back.transformed$pred,
                   ecec_sel$base_lmm_cv$back.transformed$pred,
                   ecec_sel$best_lm_cv$back.transformed$pred,
                   ecec_sel$best_lmm_cv$back.transformed$pred)
str(pred)
pev <- data.frame(clay_sel$base_lm_cv$back.transformed$pev,
                  clay_sel$base_lmm_cv$back.transformed$pev,
                  clay_sel$best_lm_cv$back.transformed$pev,
                  clay_sel$best_lmm_cv$back.transformed$pev,
                  orca_sel$base_lm_cv$back.transformed$pev,
                  orca_sel$base_lmm_cv$back.transformed$pev,
                  orca_sel$best_lm_cv$back.transformed$pev,
                  orca_sel$best_lmm_cv$back.transformed$pev,
                  ecec_sel$base_lm_cv$back.transformed$pev,
                  ecec_sel$base_lmm_cv$back.transformed$pev,
                  ecec_sel$best_lm_cv$back.transformed$pev,
                  ecec_sel$best_lmm_cv$back.transformed$pev)
str(pev)
models <- c("clay_sel$base_lm", "clay_sel$base_lmm", "clay_sel$best_lm",
            "clay_sel$best_lmm", "orca_sel$base_lm", "orca_sel$base_lmm",
            "orca_sel$best_lm", "orca_sel$best_lmm", "ecec_sel$base_lm",
            "ecec_sel$base_lmm", "ecec_sel$best_lm", "ecec_sel$best_lmm")
stats <- list()
for (i in 1:dim(obs)[2]) {
  stats[[i]] <- cvStats(obs[, i], pred[, i], pev[, i])
}
nam <- names(stats[[1]])
stats <- t(matrix(unlist(stats), ncol = 12))
colnames(stats) <- nam
stats <- data.frame(stats)
Structure <- rep(c("LRM", "LMM"), 6)
Model <- data.frame(Structure, stats[, c("me", "mae", "rmse", "srmse", "r2")])
rgroup <- c("CLAY (g kg$^{-1}$)", "SOC (g kg$^{-1}$)", "ECEC (mmol kg$^{-1}$)")
colheads <- c("Type", "ME", "MAE", "RMSE", "SRMSE", "AVE")
rowname <- rep(c("Baseline", "", "Best performing", ""), 3)
long_cap <- paste("Statistics$^a$ of the LOO-CV of baseline and best ",
                  "performing multiple linear regression models (LRMs) and ",
                  "linear mixed models (LMMs).", sep = "")
foot <- paste("${}^a$ Statistics: mean error (\\textit{ME}), mean absolute ",
              "error (\\textit{MAE}), root mean squared error ",
              "(\\textit{RMSE}), scaled root mean squared error ",
              "(\\textit{SRMSE}, unitless), and amount of variance explained ", 
              "(\\textit{AVE}, percent).", sep = "")
file <-  paste(tab_dir, "TAB05.tex", sep = "")
digits <- c(0, 2, 1, 1, 2, 1)
latex(Model, ctable = TRUE, n.rgroup = c(4, 4, 4), rgroup = rgroup, 
      file = file, label = "tab:cv-stats", table.env = TRUE, 
      cgroupTexCmd = NULL, rgroupTexCmd = NULL, rowname = rowname,
      colheads = colheads, na.blank = TRUE, caption = long_cap, where = NULL,
      size = "scriptsize", insert.bottom = foot, cdec = digits)

# SPATIAL PREDICTION - KRIGING #################################################
# Spatial predictions at unvisited locations are done using universal kriging,
# also known as kriging with external drift (trend). The entire prediction step 
# takes some time to complete because of the back-transformation of predicted
# values which is done using simulation. Because we are dealing with a large 
# amount of prediction locations, the prediciton grid has to be split into
# several tiles. Otherwise R craches due to memory overload.

# prepare base data
system("r.mask -o input=buffer_BASIN_10")
locations <- readRAST6("dnos.raster")

# CLAY - base linear mixed model -----------------------------------------------
# prepare base data
model <- clay_sel$base_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- clay_sel$base_lmm
file <- paste(r_data, "1stArticle_clay_base_krige.rda", sep = "")
covars  <- colnames(geodata$covariate)
trend   <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)

# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CLAY - best linear mixed model -----------------------------------------------
# prepare base data
model <- clay_sel$best_lmm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- clay_sel$best_lmm
file <- paste(r_data, "1stArticle_clay_best_krige.rda", sep = "")
covars <- colnames(geodata$covariate)
trend <- formula(paste("~", paste(covars, collapse = " + ")))
covars  <- readRAST6(covars)
# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars,
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CLAY - kriging predictions ---------------------------------------------------

# load base model predictions
load(paste(r_data, "1stArticle_clay_base_krige.rda", sep = ""))
clay_base_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(clay_base_lmm_krige) <- ~ x.coord + y.coord
gridded(clay_base_lmm_krige) <- TRUE
proj4string(clay_base_lmm_krige) <- wgs1984utm22s

# load best model predictions
load(paste(r_data, "1stArticle_clay_best_krige.rda", sep = ""))
clay_best_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(clay_best_lmm_krige) <- ~ x.coord + y.coord
gridded(clay_best_lmm_krige) <- TRUE
proj4string(clay_best_lmm_krige) <- wgs1984utm22s

# prepare plots with predictions
breaks <- seq(min(cal_data$CLAY), max(cal_data$CLAY), by = 1)
max_z <- max(clay_base_lmm_krige$krige.pred, clay_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
psd.colors <- colorRampPalette(R_pal$tex_pal)
col <- psd.colors(length(breaks)-1)
clay_base_lmm_krige$krige.pred.cut <- cut(clay_base_lmm_krige$krige.pred, 
                                          breaks = breaks)
clay_best_lmm_krige$krige.pred.cut <- cut(clay_best_lmm_krige$krige.pred, 
                                          breaks = breaks)
p1 <- spplot(clay_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(clay_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with predictions
dev.off()
pdf(file = paste(fig_dir, "FIG07a.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG07d.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# prepare plots with prediction variance
# take the square root to improve plotting
clay_base_lmm_krige$krige.var <- sqrt(clay_base_lmm_krige$krige.var)
clay_best_lmm_krige$krige.var <- sqrt(clay_best_lmm_krige$krige.var)
max_z <- c(max(clay_base_lmm_krige$krige.var), 
           max(clay_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(clay_base_lmm_krige$krige.var, clay_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
clay_base_lmm_krige$krige.var.cut <- cut(clay_base_lmm_krige$krige.var, 
                                         breaks = breaks)
clay_best_lmm_krige$krige.var.cut <- cut(clay_best_lmm_krige$krige.var, 
                                         breaks = breaks)
p1 <- spplot(clay_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(clay_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with prediction variance
dev.off()
pdf(file = paste(fig_dir, "FIG08a.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG08d.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# CARBON - base linear mixed model ---------------------------------------------
# prepare base data

model <- orca_sel$base_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- orca_sel$base_lmm
file <- paste(r_data, "1stArticle_orca_base_krige.rda", sep = "")
covars  <- colnames(geodata$covariate)
trend <- formula(paste("~", paste(covars, collapse = " + ")))
covars <- readRAST6(covars)

# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CARBON - best linear mixed model ---------------------------------------------
# prepare base data

model <- orca_sel$best_lm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- orca_sel$best_lmm
file <- paste(r_data, "1stArticle_orca_best_krige.rda", sep = "")
covars <- colnames(geodata$covariate)
trend <- formula(paste("~", paste(covars, collapse = " + ")))
covars <- readRAST6(covars)

# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# CARBON - kriging predictions -------------------------------------------------

# load base model predictions
load(paste(r_data, "1stArticle_orca_base_krige.rda", sep = ""))
carbon_base_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(carbon_base_lmm_krige) <- ~ x.coord + y.coord
gridded(carbon_base_lmm_krige) <- TRUE
proj4string(carbon_base_lmm_krige) <- wgs1984utm22s

# load best model predictions
load(paste(r_data, "1stArticle_orca_best_krige.rda", sep = ""))
carbon_best_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(carbon_best_lmm_krige) <- ~ x.coord + y.coord
gridded(carbon_best_lmm_krige) <- TRUE
proj4string(carbon_best_lmm_krige) <- wgs1984utm22s

# prepare plots with predictions
breaks <- seq(min(cal_data$ORCA), max(cal_data$ORCA), by = 1)
max_z <- max(carbon_base_lmm_krige$krige.pred, carbon_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
soc.colors <- colorRampPalette(R_pal$soc_pal)
col <- soc.colors(length(breaks)-1)
carbon_base_lmm_krige$krige.pred.cut <- cut(carbon_base_lmm_krige$krige.pred, 
                                            breaks = breaks)
carbon_best_lmm_krige$krige.pred.cut <- cut(carbon_best_lmm_krige$krige.pred, 
                                            breaks = breaks)
p1 <- spplot(carbon_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(carbon_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with predictions
dev.off()
pdf(file = paste(fig_dir, "FIG07b.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG07e.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# prepare plots with prediction variance
# take the square root to improve plotting
carbon_base_lmm_krige$krige.var <- sqrt(carbon_base_lmm_krige$krige.var)
carbon_best_lmm_krige$krige.var <- sqrt(carbon_best_lmm_krige$krige.var)
max_z <- c(max(carbon_base_lmm_krige$krige.var), 
           max(carbon_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(carbon_base_lmm_krige$krige.var, carbon_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
carbon_base_lmm_krige$krige.var.cut <- cut(carbon_base_lmm_krige$krige.var, 
                                           breaks = breaks)
carbon_best_lmm_krige$krige.var.cut <- cut(carbon_best_lmm_krige$krige.var, 
                                           breaks = breaks)
p1 <- spplot(carbon_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(carbon_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with prediction variance
dev.off()
pdf(file = paste(fig_dir, "FIG08b.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG08e.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# Check extrapolations in the predictions
a <- which(carbon_best_lmm_krige$krige.pred >= max(cal_data$carbon))
b <- coordinates(carbon_best_lmm_krige[a, ])
pts <- data.frame(b, carbon_best_lmm_krige$krige.pred[a])
coordinates(pts) <- ~ x.coord + y.coord
covar <- raster(readRAST6("TPI_10_15"))
image(covar, asp = 1)
points(pts, cex = 0.1, pch = 20)
pts$TPI_10_15 <- extract(x = covar, y = pts)
quantile(round(cal_data$TPI_10_15), prob = seq(0, 1, 0.1))
quantile(round(pts$TPI_10_15), prob = seq(0, 1, 0.1))

# ECEC - base linear mixed model -----------------------------------------------

# prepare base data
model <- ecec_sel$base_lmm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- ecec_sel$base_lmm
file <- paste(r_data, "1stArticle_ecec_base_krige.rda", sep = "")
covars <- colnames(geodata$covariate)
trend <- formula(paste("~", paste(covars, collapse = " + ")))
covars <- readRAST6(covars)

# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# ECEC - best linear mixed model -----------------------------------------------

# prepare base data
model <- ecec_sel$best_lmm
geodata <- as.geodata(cal_data, data.col = colnames(model$model)[1], 
                      covar.col = colnames(model$model)[-1])
model <- ecec_sel$best_lmm
file <- paste(r_data, "1stArticle_ecec_best_krige.rda", sep = "")
covars <- colnames(geodata$covariate)
trend <- formula(paste("~", paste(covars, collapse = " + ")))
covars <- readRAST6(covars)

# make predictions
spredict(model = model, geodata = geodata, file = file, covars = covars, 
         simul.back = TRUE, n.sim = 20000, pred.loc = locations,
         n.tiles = 1000, email = "alessandrosamuelrosa@gmail.com")

# ECEC - kriging predictions ---------------------------------------------------

# load base model predictions
load(paste(r_data, "1stArticle_ecec_base_krige.rda", sep = ""))
ecec_base_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(ecec_base_lmm_krige) <- ~ x.coord + y.coord
gridded(ecec_base_lmm_krige) <- TRUE
proj4string(ecec_base_lmm_krige) <- wgs1984utm22s

# load best model predictions
load("1stArticle_ecec_best_krige.rda")
ecec_best_lmm_krige <- data.frame(pred[1:4])
rm(pred)
coordinates(ecec_best_lmm_krige) <- ~ x.coord + y.coord
gridded(ecec_best_lmm_krige) <- TRUE
proj4string(ecec_best_lmm_krige) <- wgs1984utm22s

# prepare plots with predictions
breaks <- seq(min(cal_data$ECEC), max(cal_data$ECEC), by = 1)
max_z <- max(ecec_base_lmm_krige$krige.pred, ecec_best_lmm_krige$krige.pred)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
ecec.colors <- colorRampPalette(R_pal$CEC_pal)
col <- ecec.colors(length(breaks)-1)
ecec_base_lmm_krige$krige.pred.cut <- cut(ecec_base_lmm_krige$krige.pred, 
                                          breaks = breaks)
ecec_best_lmm_krige$krige.pred.cut <- cut(ecec_best_lmm_krige$krige.pred, 
                                          breaks = breaks)
p1 <- spplot(ecec_base_lmm_krige, "krige.pred.cut", col.regions = col)
p2 <- spplot(ecec_best_lmm_krige, "krige.pred.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with predictions
dev.off()
pdf(file = paste(fig_dir, "FIG07c.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG07f.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# prepare plots with prediction variance
ecec_base_lmm_krige$krige.var <- sqrt(ecec_base_lmm_krige$krige.var)
ecec_best_lmm_krige$krige.var <- sqrt(ecec_best_lmm_krige$krige.var)
max_z <- c(max(ecec_base_lmm_krige$krige.var), 
           max(ecec_best_lmm_krige$krige.var))
max_z <- max_z[which.min(max_z)]
breaks <- seq(0, max_z, by = 1)
max_z <- max(ecec_base_lmm_krige$krige.var, ecec_best_lmm_krige$krige.var)
if (max(breaks) < max_z) {
  breaks <- c(breaks, max_z)
}
col <- bpy.colors(length(breaks)-1)
ecec_base_lmm_krige$krige.var.cut <- cut(ecec_base_lmm_krige$krige.var, 
                                         breaks = breaks)
ecec_best_lmm_krige$krige.var.cut <- cut(ecec_best_lmm_krige$krige.var, 
                                         breaks = breaks)
p1 <- spplot(ecec_base_lmm_krige, "krige.var.cut", col.regions = col)
p2 <- spplot(ecec_best_lmm_krige, "krige.var.cut", col.regions = col)
p1$legend$right$args$key <- p1$legend$right$args$key[1:3]
p2$legend$right$args$key <- p2$legend$right$args$key[1:3]

# save plots with prediction variance
dev.off()
pdf(file = paste(fig_dir, "FIG08c.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p1)
dev.off()
pdf(file = paste(fig_dir, "FIG08f.pdf", sep = ""),
    width = 6.3/cm(1), height = 6.3/cm(1))
trellis.par.set(fontsize = list(text = 7, points = 5),
                plot.line = list(lwd = 0.001),
                axis.line = list(lwd = 0.01),
                layout.widths = list(left.padding = 0, right.padding = 0), 
                layout.heights = list(top.padding = 0, bottom.padding = 0))
plot(p2)
dev.off()
rm(p1, p2)
gc()

# SAVE DATA ####################################################################
ls()
# selected models and geodata
save(clay_base_geodata, clay_best_geodata, clay_base_lmm, clay_base_vario,
     clay_base_lm, clay_best_lm, clay_best_lmm, clay_best_vario,
     carbon_base_geodata, carbon_best_geodata, carbon_base_lmm, 
     carbon_base_vario, carbon_base_lm, carbon_best_lm, carbon_best_lmm,
     carbon_best_vario,
     ecec_base_geodata, ecec_best_geodata, ecec_base_lmm, ecec_base_vario, 
     ecec_base_lm, ecec_best_lm, ecec_best_lmm, ecec_best_vario,
     clay_sel, soc_sel, ecec_sel,
     file = "sm-dnos-phd-chap1-final-models.RData")

# other data
save(soil1, soil2, land1, land2, geo1, geo2, dem1, dem2, sat1, sat2,
     combs, preds, forms, drop, deltaR2,
     # base models
     clay_full, clay_vif, clay_both, clay_forward, clay_backward, 
     carbon_full, carbon_vif, carbon_both, carbon_forward, carbon_backward, 
     ecec_full, ecec_vif, ecec_both, ecec_forward, ecec_backward, 
     # base models stats
     clay_full_stats, clay_vif_stats, clay_both_stats, clay_forward_stats,
     clay_backward_stats, carbon_full_stats, carbon_vif_stats, carbon_both_stats, 
     carbon_forward_stats, carbon_backward_stats, ecec_full_stats, 
     ecec_vif_stats, ecec_both_stats, ecec_forward_stats, ecec_backward_stats,
     # cross-validation
     clay_base_lm_cv,    clay_base_lm_cv_cdf,
     clay_base_lmm_cv,   clay_base_lmm_cv_cdf,
     clay_best_lm_cv,    clay_best_lm_cv_cdf,
     clay_best_lmm_cv,   clay_best_lmm_cv_cdf,
     carbon_base_lm_cv,  carbon_base_lm_cv_cdf,
     carbon_base_lmm_cv, carbon_base_lmm_cv_cdf,
     carbon_best_lm_cv,  carbon_best_lm_cv_cdf,
     carbon_best_lmm_cv, carbon_best_lmm_cv_cdf,
     ecec_base_lm_cv,    ecec_base_lm_cv_cdf,
     ecec_base_lmm_cv,   ecec_base_lmm_cv_cdf,
     ecec_best_lm_cv,    ecec_best_lm_cv_cdf,
     ecec_best_lmm_cv,   ecec_best_lmm_cv_cdf,
     file = "sm-dnos-phd-chap1.RData")
# End!
