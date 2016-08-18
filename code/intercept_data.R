library(sp)
library(rgeos)
library(GISTools)
library(LatticeKrig)
load("../code/RData/Spatial911PtPtrn.RData")
load("../code/RData/AllTempData.RData")

FacToNum <- function(x) {
  as.numeric(as.character(x))
}

sp_pts <- SpatialPoints(cbind(pp.grid[,2], pp.grid[,1]))
sp_grd <- points2grid(sp_pts, tolerance=0.000763359)
grd_layer <- as.SpatialPolygons.GridTopology(sp_grd)
grd_cont <- ifelse(colSums(gContains(grd_layer, sp_pts, byid=T)) > 0, TRUE, FALSE)
grd_layer <- grd_layer[grd_cont]
grd_layer <- SpatialPolygonsDataFrame(grd_layer,
                                      data = data.frame(c(1:length(sp_pts))),
                                      match.ID = FALSE)
names(grd_layer) <- "ID"
row.names(grd_layer) <- paste("g", 1:length(sp_pts), sep="")

# Merging AC data with 2000 census data
AC_shp <- readShapePoly("../data/SIMMER_PARCEL_BGs/Census_2010_BG_2.shp")
census_shp <- readShapePoly("../data/Census2000/HoustonCityLimitsCensu2000.shp")
row.names(AC_shp) <- paste("AC", AC_shp$OBJECTID, sep = "")
row.names(census_shp) <- paste("CEN", census_shp$ObjectID, sep = "")

int_shp <- gIntersection(AC_shp, census_shp, byid = T, drop_lower_td = T)
int_id <- strsplit(names(int_shp), " ")

AC_id <- (sapply(int_id, "[[", 1))
census_id <- (sapply(int_id, "[[", 2))

int_area <- gArea(int_shp, byid = T)
AC_area <- gArea(AC_shp, byid = T)
census_area <- gArea(census_shp, byid = T)

# This method, where you use the proportion of the census tract area, only works for count data, not percentages
AC_index <- match(AC_id, row.names(AC_shp))
AC_areas <- AC_area[AC_index]
AC_prop <- zapsmall(int_area / AC_areas, 3)
AC_int <- zapsmall(AC_shp$NOAC[AC_index] * AC_prop, 5)



# Use different denominator to realign percentage
AC_shp$pctNOAC <- AC_shp$NOAC / AC_shp$Parcels
census_index <- match(census_id, row.names(census_shp))
census_areas <- census_area[census_index]
census_prop <- zapsmall(int_area / census_areas)
pctAC_int <- zapsmall(AC_shp$pctNOAC[AC_index] * census_prop, 5)
df <- data.frame(census_id, AC_id, pctAC_int, AC_int)

pctNOAC_realign <- xtabs(df$pctAC_int ~ df$census_id)
NOAC_realign <- xtabs(df$AC_int ~ df$census_id)

cens_ind <- match(names(NOAC_realign), paste("CEN", census_shp$ObjectID, sep = ""))
census_shp$NOAC <- numeric(nrow(census_shp))
census_shp$NOAC[1:nrow(census_shp)] <- NA
census_shp$NOAC[cens_ind] <- NOAC_realign
# census_shp$NOAC[is.na(census_shp$NOAC)] <- -1
census_shp$pctNOAC <- numeric(nrow(census_shp))
census_shp$pctNOAC[1:nrow(census_shp)] <- NA
census_shp$pctNOAC[cens_ind] <- pctNOAC_realign
# census_shp$pctNOAC[is.na(census_shp$pctNOAC)] <- -1

# AC_cont <- ifelse(colSums(gIntersects(AC_shp, grd_layer, byid=T)) > 0, TRUE, FALSE)
# AC_shp2 <- AC_shp[AC_cont, ]

# Eye test to make sure things make sense
# par(mfrow = c(1, 2))
# par(mar=c(0, 0, 0, 0))
# shades <- auto.shading(census_shp$NOAC, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(census_shp, census_shp$NOAC, shades)
# choro.legend(-95.7, 30.3, shades, title="NOAC", cex=0.8)
# shades <- auto.shading(AC_shp2$NOAC, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(AC_shp2, AC_shp2$NOAC, shades)
# choro.legend(-95.7, 30.275, shades, title="NOAC", cex=0.8)
# 
# par(mfrow = c(1, 2))
# par(mar=c(0, 0, 0, 0))
# shades <- auto.shading(census_shp$pctNOAC, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(census_shp, census_shp$pctNOAC, shades)
# choro.legend(-95.7, 30.3, shades, title="%NOAC", cex=0.8)
# shades <- auto.shading(AC_shp2$pctNOAC, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(AC_shp2, AC_shp2$pctNOAC, shades)
# choro.legend(-95.7, 30.275, shades, title="%NOAC", cex=0.8)

###############################################################################
### Align census data with grid                                             ###
###############################################################################

grid_int <- gIntersection(grd_layer, census_shp, byid=T)

tmp <- strsplit(names(grid_int), " ")
grid_id <- (sapply(tmp, "[[", 1))
census_id <- (sapply(tmp, "[[", 2))

int_area <- gArea(grid_int, byid=T)
cen_area <- gArea(census_shp, byid=T)
grid_area <- gArea(grd_layer, byid=T)

# This method, where you use the proportion of the census tract area, only works for count data, not percentages
c_ind <- match(census_id, row.names(census_shp))
census_areas <- cen_area[c_ind]
census_prop <- zapsmall(int_area/census_areas, 3)
pop <- zapsmall(census_shp$Total[c_ind] * census_prop, 5)
NOAC <- zapsmall(census_shp$NOAC[c_ind] * census_prop, 5)

# Align percentages to grid
grid_areas <- grid_area[grid_id]
grid_prop <- zapsmall(int_area / grid_areas, 3)
over65PCT <- zapsmall(census_shp$over65PCT[c_ind] * grid_prop, 5)
NOACPCT <- zapsmall(census_shp$pctNOAC[c_ind] * grid_prop, 5)
disabledPCT <- zapsmall(census_shp$disabledPC[c_ind] * grid_prop, 5)
HispanicPCT <- zapsmall(census_shp$HispanicPC[c_ind] * grid_prop, 5)
BlackPCT <- zapsmall(census_shp$BlackPCT[c_ind] * grid_prop, 5)
under5PCT <- zapsmall(census_shp$under5PCT[c_ind] * grid_prop, 5)
alonePCT <- zapsmall(census_shp$alonePCT[c_ind] * grid_prop, 5)
povertyPCT <- zapsmall(census_shp$povertyPCT[c_ind] * grid_prop, 5)


df <- data.frame(grid_id, census_id, TotalPopulation = pop, over65PCT,
                 under5PCT, BlackPCT, HispanicPCT, disabledPCT, 
                 alonePCT, povertyPCT, NOACPCT, NOAC)

int_over65PCT <- xtabs(df$over65PCT ~ df$grid_id)
int_under5PCT <- xtabs(df$under5PCT ~ df$grid_id)
int_pop <- xtabs(df$TotalPopulation ~ df$grid_id)
int_NOAC <- xtabs(df$NOAC ~ df$grid_id)
int_NOACPCT <- xtabs(df$NOACPCT ~ df$grid_id)
int_BlackPCT <- xtabs(df$BlackPCT ~ df$grid_id)
int_HispanicPCT <- xtabs(df$HispanicPCT ~ df$grid_id)
int_disabledPCT <- xtabs(df$disabledPCT ~ df$grid_id)
int_alonePCT <- xtabs(df$alonePCT ~ df$grid_id)
int_povertyPCT <- xtabs(df$povertyPCT ~ df$grid_id)

index <- as.numeric(gsub("g", "", names(int_over65PCT)))
TotalPopulation <- over65PCT <- under5PCT <- NOAC <- NOACPCT <- BlackPCT <-
  HispanicPCT <- disabledPCT <- alonePCT <- povertyPCT <- 
  as.numeric(rep(NA, nrow(grd_layer)))

TotalPopulation[index] <- int_pop
over65PCT[index] <- int_over65PCT
under5PCT[index] <- int_under5PCT
NOAC[index] <- int_NOAC
NOACPCT[index] <- int_NOACPCT
BlackPCT[index] <- int_BlackPCT
HispanicPCT[index] <- int_HispanicPCT
disabledPCT[index] <- int_disabledPCT
alonePCT[index] <- int_alonePCT
povertyPCT[index] <- int_povertyPCT
data_grid <- SpatialPolygonsDataFrame(grd_layer,
                data = data.frame(data.frame(grd_layer), TotalPopulation, 
                                  over65PCT, under5PCT, NOAC, NOACPCT,
                                  BlackPCT, HispanicPCT, disabledPCT, alonePCT,
                                  povertyPCT, coordinates(grd_layer)), 
                match.ID = F)




pct65[index] <- int.layer.pct65
pop[index] <- int.layer.pop
noac[index] <- int.layer.noac
pct_noac[index] <- int.layer.pct_noac
pct_under5[index] <- int.layer.pct_under5
pct_black[index] <- int.layer.pct_black
pct_hispanic[index] <- int.layer.pct_hispanic
sql[index] <- int.layer.sql
int.layer <- SpatialPolygonsDataFrame(grd.layer, 
                                      data=data.frame(data.frame(grd.layer), 
                                                      pop, pct65, pct_under5, noac, 
                                                      pct_noac, pct_black, pct_hispanic,
                                                      sql, coordinates(grd.layer)), match.ID=FALSE)

# png("TotalPopulation.png", width=720, height=540)
par(mfrow = c(1, 2))
par(mar=c(0, 0, 0, 0))
shades <- auto.shading(census_shp$Total, n = 9, cols=brewer.pal(9, "Greens"))
choropleth(census_shp, census_shp$Total, shades)
choro.legend(-95.7, 30.3, shades, title="Total Population", cex=0.8)
shades <- auto.shading(data_grid$TotalPopulation, n = 9, cols=brewer.pal(9, "Greens"))
choropleth(data_grid, data_grid$TotalPopulation, shades)
choro.legend(-95.7, 30.275, shades, title="Total Population", cex=0.8)
# dev.off()

# png("PercentOver65.png", width=720, height=540)
par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0))
cuts <- quantile(census_shp$pctNOAC, probs = seq(0, 1, length=10), na.rm = T)
cuts <- cuts[-c(1, 10)]
shades <- auto.shading(census_shp$pctNOAC, cutter = function(x, n, params) {cuts}, n = 9, cols=brewer.pal(9, "Reds"))
choropleth(census_shp, census_shp$pctNOAC, shades)
choro.legend(-95.7, 30.3, shades, title="Percent Over 65", cex=0.8)
shades <- auto.shading(data_grid$NOACPCT, n = 9, cols = brewer.pal(9, "Reds"))
choropleth(data_grid, data_grid$NOACPCT, shades)
choro.legend(-95.7, 30.275, shades, title="Percent Over 65", cex=0.8)
# dev.off()

# Convert into a useful data.frame
grid_df <- data.frame(data_grid)
names(grid_df)[c(12,13)]  <- c("Longitude", "Latitude")
grid_df <- grid_df[order(grid_df$Longitude, grid_df$Latitude), ]

vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
getSpatialMeans <- function(temp.var) {
  H <- lapply(temp.data.nomiss, function(x) { as.matrix(x[, temp.var]) })
  spatial.means <- apply(sapply(H, function(x) { apply(x, 1, mean) }), 1, mean)
  spatial.means
}

spatial_means <- lapply(vars, getSpatialMeans)
names(spatial_means) <- vars

intercept_df <- cbind(pp.grid, 
                      grid_df[, !(names(grid_df) %in% c("Longitude", "Latitude"))], 
                      spatial_means)
save(intercept_df, file="./RData/RevisedInterceptData.RData")


