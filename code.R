## Packages

## The downscaling described in this paper has been implemented using
## the free software environment R and several contributed packages:

## - =raster= for spatial data manipulation and analysis.

## - =solaR= for the solar geometry.

## - =gstat= and =sp= for the geostatistical analysis.

## - =parallel= for multi-core parallelization.

## - =rasterVis= for spatial data visualization methods.


library(sp)
library(raster)
library(maptools)
library(gstat)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(solaR)
library(parallel)

## Change 'MY_FOLDER' with the name of the folder
## in which the repository has been cloned.
setwd('MY_FOLDER')

##################################################################
## Radiation Data
##################################################################

## Satellite data can be freely downloaded from CM SAF (www.cmsaf.eu),
## previous registration and log-in, choosing hourly climate data sets
## named /SIS/ (Global Horizontal Irradiation)) and /SID/ (Beam
## Horizontal Irradiation) for 2005. Due to the size of these files,
## they have not been uploaded it into this github repository.
## Instead, the resulting SISa2005 raster is uploaded in the github
## repository. 
SISa2005 <- raster('data/SISa2005')
## Beware that the projection of SISa2005 is long-lat. This projection is needed to compute the sun geometry. After that step, this raster will be projected to the UTM projection.

## DO NOT RUN unless you have downloaded the hourly CMSAF files
## Once downloaded the irradiation files, they are collated into a single RasterStack.
## GHI
listFich <- dir(pattern='SIShm2005')
stackSIS <- stack(listFich)
## BHI
listFich <- dir(pattern='SIDhm2005')
stackSID <- stack(listFich)

## Annual GHI: it is computed *before* the UTM projection
SISa2005 <- calc(stackSIS,sum,na.rm=TRUE)

## Both rasters are projected to the UTM projection for compatibility with the Digital Elevation Model. 
projUTM  <-  CRS('+proj=utm +zone=30')

stackSIS <- projectRaster(stackSIS, crs=projUTM)
stackSID <- projectRaster(stackSID, crs=projUTM)
## END of DO NOT RUN

##################################################################
## Digital Elevation Model
##################################################################

## The DEM is obtained from www.ign.es MDT-200 files.  This DEM uses
## the UTM projection. Only the DEM corresponding to La Rioja area is
## needed:

elevRioja <- raster('data/elev')
names(elevRioja) <- 'elev'

##################################################################
## Sun geometry
##################################################################

## The first step is to compute the sun angles (height and azimuth)
## and the extraterrestial solar irradiation for every cell of the
## CMSAF rasters. The function =calcSol= from the =solaR= package=
## calculates the sun angles. Its results, combined with the
## =overlay= method from the =raster= package, produce three
## multilayer =Raster= objects with the sun geometry needed for the
## next steps. 

## The calculation of sun geometry is performed with the resolution of
## CM SAF.
latlon <- stack(init(SISa2005, v='y'),
                init(SISa2005, v='x'))
names(latlon) <- c('lat', 'lon')

## Sun angles are calculated with 5 min samples
BTi <- seq(as.POSIXct('2005-01-01 00:00:00'),
           as.POSIXct('2005-12-31 23:55:00'), by='5 min')
## But are aggregated to hourly samples with this function, that
## extract *and center* hour for aggregation
hour <- function(tt)as.POSIXct(trunc(tt + 30*60, 'hours'))
## Thus, this is the hourly time index
BTh <- unique(hour(BTi))
############################################################
## Extraterrestial solar irradiation
############################################################

## The SIS coordinates (locs) are traversed with mclapply so for
## every point of the Raster object a time series of sun geometry is
## computed
Bo05min <- overlay(latlon, fun=function(lat, lon)
{
    locs <- as.data.frame(rbind(lat, lon))
    b <- mclapply(locs, function(p)
    {
        cat('Point: ', p[2], ' ', p[1], '\n')
        ## Mean solar Time
        hh <- local2Solar(BTi, p[2])
        ## Sun geometry
        sol <- calcSol(p[1], BTi=hh)
        ## Extraterrestial solar irradiation
        as.data.frameI(sol)$Bo0
    }, mc.cores = detectCores())
    res <- do.call(rbind, b)
})
## Bo05min is RasterBrick object with a layer for each element of
## the time index BTi. It can be set with setZ
Bo05min <- setZ(Bo05min, BTi)
names(Bo05min) <- as.character(BTi)

## Both zApply and projectRaster accept parallel computing
beginCluster(type = 'PSOCK')
## Hourly aggregation with zApply
Bo0h <- clusterR(Bo05min, zApply, args = list(by = hour))
## clusterR does not preserve the slot z set by zApply
Bo0h <- setZ(Bo0h, BTh)
## Project to UTM (same as DEM raster)
Bo0h <- projectRaster(Bo0h, crs=projUTM,
                     filename = 'Bo0h', overwrite = TRUE)
## We don't need the cluster anymore
endCluster()

############################################################
## Sun height
############################################################

AlS <- overlay(latlon, fun=function(lat, lon)
{
  locs <- as.data.frame(rbind(lat, lon))
  b <- mclapply(locs, function(p)
  {
      cat('Point: ', p[2], ' ', p[1], '\n')
      hh <- local2Solar(BTi, p[2])
      sol <- calcSol(p[1], BTi=hh)
      ## Calculation of the sun height and rounding with 2
      ## decimal positions to reduce the matrix size.
      round(r2d(as.data.frameI(sol)$AlS), 2)
  }, mc.cores = detectCores())
  res <- do.call(rbind, b)
})

AlS <- setZ(AlS, BTi)
names(AlS) <- as.character(BTi)

beginCluster(n = detectCores(), type = 'PSOCK')
AlSh <- clusterR(AlS, zApply, args = list(by = hour))
AlSh <- setZ(AlSh, BTh)
AlSh <- projectRaster(AlSh, crs=projUTM,
                      filename = 'AlSh', overwrite = TRUE)
endCluster()

############################################################
## Azimuth
############################################################

AzS <- overlay(latlon, fun=function(lat, lon){
  locs <- as.data.frame(rbind(lat, lon))
  b <- mclapply(locs, function(p)
  {
      cat('Point: ', p[2], ' ', p[1], '\n')
      hh <- local2Solar(BTi, p[2])
      ## Calculation of the solar azimuth and rounding with 2 decimal
      ## positions to reduce the matrix size.
      sol <- calcSol(p[1], BTi=hh)
      round(r2d(as.data.frameI(sol)$AzS),2)
  }, mc.cores = detectCores())
  res <- do.call(rbind, b)
})

AzS <- setZ(AzS, BTi)
names(AzS) <- as.character(BTi)

beginCluster(n = detectCores(), type = 'PSOCK')
AzSh <- clusterR(AzS, zApply, args = list(by = hour))
AzSh <- setZ(AzSh, BTh)
AzSh <- projectRaster(AzSh, crs=projUTM,
                      filename = 'AzSh', overwrite = TRUE)
endCluster()

##################################################################
## Irradiation Components
##################################################################

## The CMSAF rasters must be transformed to the higher resolution of
## the DEM (UTM 200mx200m). As a consequence of the different pixel
## geometry between DEM (square) and irradiation rasters (rectangle)
## the process is performed in two steps. The first step increases
## the spatial resolution of the irradiation rasters to a similar and
## also larger pixel size than the DEM with =disaggregate=. The
## second step post-processes the previous step by means of a
## bilineal interpolation which resamples the raster layer and
## achieves same DEM resolution.

## On the other hand, the diffuse irradiation is obtained from the
## global and beam irradiation rasters. The two components of the
## diffuse irradiation, isotropic and anisotropic, can be separated
## with the anisotropy index, computed as the ratio between beam
## and extraterrestial irradiation.

## Scale factor
sf <- res(stackSID)/res(elev)

## extraterrestrial irradiation and sun angles
Bo0hd <- disaggregate(Bo0h, sf)
Bo0hdr <- resample(Bo0hd, elev)

AzShd <- disaggregate(AzSh, sf)
AzShdr <- resample(AzShdr, elev)

AlShd <- disaggregate(AlSh, sf)
AlShdr <- resample(AlShd, elev)

## Beam irradiation
SIDd <- disaggregate(stackSID, sf)
SIDdr <- resample(SIDd, elev)

## Global irradiation
SISd <- disaggregate(stackSIS, sf)
SISdr <- resample(SISd, elev)

## Diffuse irradiation
## k1 is the anisotropy index
k1 <- SIDdr/Bo0hdr
Difdr <- SISdr - SIDdr
Difiso <- (1 - k1) * Difdr
Difani <- k1 * Difdr

##################################################################
## Horizon Angle
##################################################################

## The maximum horizon angle required for the horizon blocking
## analysis and also to derive the SVF is obtained with the next code.
## For each direction angle (value of the =alfa= vector) the maximum
## horizon angle is calculated for a set of points across that
## direction from each of the locations defined in =xyelev= (derived
## from the DEM raster and transformed in the matrix =locs= visited by
## rows).  Separations between the origin locations and points along
## each direction are defined in the vector =seps=. The elevation
## (=z1=) of these points is converted into the horizon angle: the
## largest of these angles is the horizon angle for that
## direction. The result of each =apply= step is a matrix which is
## used to fill in a RasterLayer (=r=). The result of =mclapply= is a
## list, =hor=, of =RasterLayer= which can be converted into a
## =RasterStack= with =stack=. Each layer of this =RasterStack=
## corresponds to a different direction.

## Maximum distance
d <- 20000

## The DEM used for the horizon computation must be smaller according
## to this distance
extHorizon <- extent(elevRioja) + c(d, -d, d, -d)
elev <- crop(elevRioja, extHorizon)

## Sample distance defined by DEM resolution
resD <- max(res(elev))
## Separations
seps <- seq(resD, d, by=resD)

## Raster definition with UTMX, UTMY coordinates and elevation
xyelev <- stack(init(elev, v='x'),
                init(elev, v='y'),
                elev)
names(xyelev) <- c('x', 'y', 'elev')
## Coordinates
locs <- as.matrix(xyelev)
px <- locs[, 1]
py <- locs[, 2]
pz <- locs[, 3]

## Angle of sector sampling 5º.
inc <- pi/36
alfa <- seq(-pi, pi - inc, inc)

## Compute the horizon angle for each direction (alfa) with parallel
## computing using mclapply.
hor <- mclapply(alfa, function(ang)
{
    cat('Angle: ', ang, '\n')
    ## Loop across separations
    hor <- sapply(seps, function(sep)
    {
        x1 <- px - sin(ang) * sep
        y1 <- py - cos(ang) * sep
        p1 <- cbind(x1,y1)
        z1 <- elevRioja[cellFromXY(elevRioja, p1)]
        r2d(atan2(z1 - pz, sep))
    })
    ## The result is a matrix of angles with a row for each point of
    ## the DEM (row of locs) and with a column for each separation.
    r <- raster(elev)
    ## The horizon angle at each location is the maximum angle, that
    ## must be positive: "apply" computes the maximum value per row. 
    r[] <- apply(hor, 1, max, 0)
    r
}, mc.cores = detectCores())
## The result is a list of RasterLayer, one per each direction angle.
## Create a RasterStack with a layer for each angle
horizon <- stack(hor)
writeRaster(horizon, 'horizon')

##################################################################
## Horizon Blocking
##################################################################

## The horizon blocking is analyzed evaluating the solar geometry in
## 15 minutes samples, particularly the solar elevation and azimuth
## angles throughout the original irradiation raster. Secondly, the
## hourly averages ared calculated, disaggregated and post-processed
## as previously explained for the irradiation rasters. The decision
## of solving the solar geometry with low resolution rasters allows
## for the significan reduction of computational time without
## penalizing results.

## Cut the Azimut raster in classes according to the alfa vector
## (directions). 
idxAngle <- cut(AzShdr, breaks=r2d(alfa))
## With idxAngle the values of the horizon Raster corresponding to
## each angle class can be extracted using stackSelect. 
AngAlt <- stackSelect(horizon, idxAngle)
## The number of layers of AngAlt is the same as idxAngle and,
## therefore, the same as AzShdr. It can be used for comparison
## with with solar height, AlS. If AngAlt is greater, there is
## horizon blocking.
dilogical <- (AngAlt - AlShdr) < 0

## Beam irradiation corrected with horizon blocking.
Dirh <- SIDdr * dilogical
## Diffuse anisotrophic irradiation corrected with horizon
## blocking
Difani <- Difani * dilogical
##################################################################
## Sky View Factor
##################################################################

## The Sky View Factor can be easily computed from the =horizon= object
## with the equation proposed above.  

## Calculation with  Ruiz-Arias et al. (2010) equation.
SVFRuizArias <- calc(horizon, function(x) sin(d2r(x))^2)
SVF <- 1 - mean(SVFRuizArias)
Difiso <- Difiso * SVF

## Global irradiation
GHIh <- Difani + Difiso + Dirh
## Annual sum of hourly global irradiation
GHI2005a <- calc(GHIh, fun=sum)

## GHI2005a stands for the downscaled irradiation raster without ked. 
## It has been uploaded into this github repository

GHI2005a <- raster('data/GHI2005a')

##################################################################
## Kriging with external drift
##################################################################

## The downscaled irradiation rasters can be improved using kriging
## with external drift. Irradiation data from on-ground
## meteorological stations is interpolated with the downscaled
## irradiation raster as explanatory variable.

load('data/Stations.RData')
UTM <- SpatialPointsDataFrame(Stations[,c(2,3)], Stations[,-c(2,3)],
                              proj4string=projUTM)


vgmCMSAF <- variogram(GHImed ~ GHIcmsaf, UTM)
fitvgmCMSAF <- fit.variogram(vgmCMSAF, vgm(model='Nug'))

gModel <- gstat(NULL, id='G0yKrig',
                formula= GHImed ~ GHIcmsaf,
                locations=UTM, model=fitvgmCMSAF)

names(GHI2005a) <- 'GHIcmsaf'
G0yKrig <- interpolate(GHI2005a, gModel, xyOnly=FALSE)

##################################################################
## Brief analysis of the results
##################################################################

## evaluation of downscaling + KED
G0yKrig_v <- extract(G0yKrig, UTM)
## evaluation of GHIcmsaf (downscaling without kriging)
GHI2005a_v <- extract(GHI2005a, UTM)
## evaluation of GHI by CM SAF
SISa2005_v <- extract(SISa2005, UTM)
## collate results.
results <- cbind(SISa2005_v, G0yKrig_v, GHI2005a_v)

MAE <- colMeans(abs(results - UTM$GHImed))
names(MAE) <- c('MAEcmsaf', 'MAEkrig', 'MAEcmsafdown')

## zonal estatistics sd
cells <- raster(SISa2005dr)

xy <- as(SISa2005dr, 'SpatialPoints')
cells[] <- cellFromXY(SISa2005, xy)

SD <- SD2 <- raster(SISa2005)
SD[] <- zonal(G0yKrig, cells, sd)[,2]
SD2[] <- zonal(GHI2005a, cells, sd)[,2]
