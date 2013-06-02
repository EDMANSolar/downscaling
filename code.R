
## Packages
## #+LABEL: sec:packages
## The downscaling described in this paper has been implemented using
## the free software environment R \citep{RDevelopmentCoreTeam2012}
## and several contributed packages: 

## - =raster= \citep{Hijmans.Etten2012} for spatial data manipulation
##   and analysis.
## - =solaR= \citep{Perpinan2012} for the solar
##   geometry.
## - =gstat= \citep{Pebesma2004} and =sp=
##   \citep{Pebesma.Bivand2005} for the geostatistical analysis.
## - =parallel= for multi-core parallelization.
## - =rasterVis= \citep{Perpinan.Hijmans2012} for spatial data
##   visualization methods.


library(sp)
library(raster)
rasterOptions(todisk=FALSE)
rasterOptions(chunksize = 1e+06, maxmemory = 1e+07)
library(maptools)
library(gstat)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(solaR)
library(parallel)

## Data
## #+LABEL: sec:data
## Satellite data can be freely downloaded from CM SAF[fn:1] choosing hourly
## climate data sets named /SIS/ (Global Horizontal Irradiation)) and
## /SID/ (Beam Horizontal Irradiation) for 2005. Both rasters are
## projected to the UTM projection for compatibility with the Digital
## Elevation Model.

projUTM  <-  CRS('+proj=utm +zone=30')
projLonLat <- CRS(' +proj=longlat +ellps=WGS84')

listFich <- dir(pattern='SIShm2005')
stackSIS <- stack(listFich)
stackSIS <- projectRaster(stackSIS,crs=projUTM)

## Annual GHI
SISa2005 <- calc(stackSIS,sum,na.rm=TRUE)

## BHI
listFich <- dir(pattern='SIDhm2005')
stackSID <- stack(listFich)
stackSID <- projectRaster(stackSID, crs=projUTM)

## The DEM provided in elevG can be crop to the region
## analyzed (La Rioja). As stated above, this DEM uses the UTM
## projection.

elevSpain <- raster('elevSpain.grd')
elev <- crop(elevSpain, extent(479600, 616200, 4639600, 4728400))
names(elev)<-'elev'

## Sun geometry
## #+LABEL: sec:sun-geometry
## The first step is to compute the sun angles (height and azimuth)
## and the extraterrestial solar irradiation for every cell of the
## CMSAF rasters. The functions =fSolD= and =fSolI= from the =solaR=
## package= calculate respectively the daily and intradaily sun
## geometry. These functions, combined with the =overlay= method from
## the =raster= package, produce three multilayer =Raster= objects
## with the sun geometry needed for the next steps. For the sake of
## brevity we show only the procedure for the extraterrestial solar
## irradiation.

## function to extract hour for aggregation
hour <- function(tt)as.POSIXct(trunc(tt, 'hours'))

## work with the resolution of CM SAF
r <- SISa2005

## Definition of a coordinate raster.
latlon <- stack(init(r, v='y'), init(r, v='x'))
names(latlon) <- c('lat', 'lon')

## Overlay permits using several layers from a RasterStack or
## RasterBrick. The calculation of sun geometry is performed with
## the resolution of CM SAF.

############################################################
## Extraterrestial solar irradiation
############################################################

## The extraterrestrial irradiation is calculated with 5 min
## samples
BTi <- seq(as.POSIXct('2005-01-01 00:00:00'),
           as.POSIXct('2005-12-31 23:55:00'), by='5 min')

B05min <- overlay(latlon, fun=function(lat, lon){
    ## every point is a column of a data.frame...
    locs <- as.data.frame(rbind(lat, lon))
    ## These columns are traversed with lapply so for every point
    ## of the Raster object a time series of sun geometry is
    ## computed
    b <- lapply(locs, function(p){
        ## Mean solar Time
        hh <- local2Solar(BTi, p[2])
        ## Sun geometry
        sol <- calcSol(p[1], BTi=hh)
        ## Extraterrestial solar irradiation
        Bo0 <- as.data.frameI(sol)$Bo0
        Bo0 })
    res <- do.call(rbind, b)})
## B05min is RasterBrick object with a layer for each element of
## the time index BTi. It can be set with setZ
B05min <- setZ(B05min, BTi)
names(B05min) <- as.character(BTi)
## Hourly aggregation with zApply
B0h <- zApply(B05min, by=hour, fun=mean)
projectRaster(B0h,crs=projUTM)

############################################################
## Sun height
############################################################

BTi <- seq(as.POSIXct('2005-01-01 00:00:00'),
           as.POSIXct('2005-12-31 23:45:00'), by='15 min')

AlS <- overlay(latlon, fun=function(lat, lon){
    locs <- as.data.frame(rbind(lat, lon))
    foo <- lapply(locs, function(p){
        ## Calculation of the local hour with local2Solar.
        hh <- local2Solar(BTi, p[2])
        sol <- calcSol(p[1], BTi=hh)
        ## Calculation of the sun height and rounding with 2
        ## decimal positions to reduce the matrix size.
        AlS <- round(r2d(as.data.frameI(sol)$AlS, 2)
        AlS})
    res <- do.call(rbind, foo)})
AlSn <- setZ(AlS, BTi)
names(AlSn) <- as.character(BTi)

AlSh <- zApply(AlSn, by=hour, fun=mean)

############################################################
## Azimuth
############################################################

BTi <- seq(as.POSIXct('2005-01-01 00:00:00'),
           as.POSIXct('2005-12-31 23:45:00'), by='15 min')

AzS <- overlay(latlon, fun=function(lat, lon){
  locs <- as.data.frame(rbind(lat, lon))
  foo <- lapply(locs, function(p){
    ## Calculation of the local hour with local2Solar.
    hh <- local2Solar(BTi, p[2])
    ## Calculation of the solar azimuth and rounding with 2 decimal
    ## positions to reduce the matrix size.
    sol <- calcSol(p[1], BTi=hh)
    AzS <- round(r2d(as.data.frameI(sol)$AzS,2)
    AzS})
  res <- do.call(rbind, foo)})
## Setting of the temporal index to AzS (15-min)
AzSn <- setZ(AzS, BTi)
names(AzSn) <- as.character(BTi)

AzSh <- zApply(AzSn, by=hour, fun=mean)

## Irradiation Components
## #+LABEL: sec:irradiation
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

## Beam irradiation
SIDd <- disaggregate(stackSID, sf)
SIDdr <- resample(SIDd, elev)

## Global irradiation
SISd <- disaggregate(stackSIS, sf)
SISdr <- resample(SISd, elev)

## Diffuse irradiation
Difdr <- SISdr-SIDdr

## extraterrestrial irradiation
B0hd <- disaggregate(B0h, sf)
B0hdr <- resample(B0hd, elev)

## anisotropy index
k1 <- SIDdr/B0hdr

Difiso <-(1-k1) * Difdr
Difani <- k1 * Difdr

## Horizon Angle
## #+LABEL: sec:horizon
## The maximum horizon angle required for the horizon blocking
## analysis and also to derive the SVF is obtained with the next
## code.  For each direction angle (value of the =alfa= vector) the
## maximum horizon angle is calculated for a set of points across
## that direction from each of the locations defined in =xyelev=
## (derived from the DEM raster and transformed in the matrix =locs=
## visited by rows).  Separations between the origin locations and points
## along each direction are defined in the vector =seps=. The elevation
## (=z1=) of these points is converted into the horizon angle: the
## largest of these angles is the horizon angle for that
## direction. The result of each =apply= step is a matrix which is used to
## fill in a RasterLayer (=r=). The result of =mclapply= is a list, =hor=,
## of =RasterLayer= which can be converted into a =RasterStack= with
## =stack=. Each layer of this =RasterStack= corresponds to a different
## direction.

## Sample distance defined by DEM resolution
resD <- max(res(elev))

## Maximum sample distance
d <- 20000
seps <- seq(resD, d, by=resD)

## Raster definition with UTMX, UTMY coordinates and elevation
xyelev <- stack(init(elev, v='x'),
                init(elev, v='y'),
                elev)
names(xyelev) <- c('x', 'y','elev')


## Angle of sector sampling 5º.
inc <- pi/36
alfa <- seq(-0.5*pi,(1.5*pi-inc), inc)

## Calculation with parallel computing using mclapply and 8 nodes.
locs <- as.matrix(xyelev)

hor <- mclapply(alfa, function(ang){
    h <- apply(locs, 1, function(p){
        x1 <- p[1]+cos(ang)*seps
        y1 <- p[2]+sin(ang)*seps
        p1 <- cbind(x1,y1)
        z1 <- elevSpain[cellFromXY(elevSpain,p1)]
        hor <- r2d(atan2(z1-p[3], seps))
        maxHor <- max(hor[which.max(hor)], 0)
    })
    r <- raster(elev)
    r[] <- matrix(h, nrow=nrow(r), byrow=TRUE)
    r}, mc.cores=8)

horizon <- stack(hor)

## Horizon Blocking
## #+LABEL: sec:block
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
idxAngle <- cut(AzShr, breaks=r2d(alfa))
## With idxAngle the values of the horizon Raster corresponding to
## each angle class can be extracted using stackSelect. 
AngAlt <- stackSelect(horizon, idxAngle)
## The number of layers of AngAlt is the same as idxAngle and,
## therefore, the same as AzShr. It can be used for comparison
## with with solar height, AlS. If AngAlt is greater, there is
## horizon blocking.
dilogical <- ((AngAlt-AlShr)<0)

## Beam irradiation corrected with horizon blocking.
Dirh <- SIDdr * dilogical
## Diffuse anisotrophic irradiation corrected with horizon
## blocking
Difani <- Difani * dilogical

## Sky View Factor
## #+LABEL: sec:svf
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

## Kriging with external drift
## #+LABEL: sec:ked
## The downscaled irradiation rasters can be improved using kriging
## with external drift. Irradiation data from on-ground
## meteorological stations is interpolated with the downscaled
## irradiation raster as explanatory variable. To define the
## variogram here we use the results previously published in
## \cite{Antonanzas-Torres.Canizares.ea2013}.

load('Stations.RData')
UTM <- SpatialPointsDataFrame(Stations[,c(2,3)], Stations[,-c(2,3)],
                              proj4string=CRS('+proj=utm +zone=30 +ellps=WGS84'))


vgmCMSAF <- variogram(GHImed~GHIcmsaf, UTM)
fitvgmCMSAF <- fit.variogram(vgmCMSAF, vgm(model='Nug'))

gModel <- gstat(NULL, id='G0yKrig',
                formula= GHImed ~ GHIcmsaf,
                locations=UTM, model=fitvgmCMSAF)

names(GHI2005a) <- 'GHIcmsaf'
G0yKrig <- interpolate(GHI2005a, gModel, xyOnly=FALSE)

## Brief analysis of the results
## #+LABEL: sec:results
## The mean absolute error ($MAE$) is analyzed for the set of
## stations (Figure~\ref{fig:mapstations}). Table~\ref{tab:MAE} shows
## the improvement in errors when considering the KED. The
## downscaling using KED presents lower MAE than with the original CM
## SAF data.

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

## #+CAPTION: MAE of different scenarios against measured data from on-ground pyranometers ($kWh/m^2$)
## #+LABEL: tab:MAE
## | $MAE_{cmsaf}$ | $MAE_{down}$ | $MAE_{KED}$ |
## |---------------+--------------+-------------|
## |        119.6 |       175.6 |       102.9 |

## On the other hand, variability due to the downscaling procedure
## can be considered with the =zonal= function from =raster= package.

## zonal estatistics sd
cells <- raster(SISa2005dr)

xy <- as(SISa2005dr, 'SpatialPoints')
cells[] <- cellFromXY(SISa2005, xy)

SD <- SD2 <- raster(SISa2005)
SD[] <- zonal(G0yKrig, cells, sd)[,2]
SD2[] <- zonal(GHI2005a, cells, sd)[,2]
