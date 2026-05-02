library(ncdf4)
fai_max<-nc_open("dataset-sargassum-cls-merged-msi-oli-global-lr.nc")
print(fai_max)


lon <- ncvar_get(fai_max,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(fai_max,"latitude")
nlat <- dim(lat)
head(lat)

time <- ncvar_get(fai_max,"time")
tunits <- ncatt_get(fai_max,"time","units")
nt <- dim(time)

fai_array<-ncvar_get(fai_max,"nfai_max")
longname<- ncatt_get(fai_max, "nfai_max", "long_name")
dunits<-ncatt_get(fai_max, "nfai_max", "units")
fillvalue<-ncatt_get(fai_max, "nfai_max", "FillValue")

title <- ncatt_get(fai_max,0,"title")
institution <- ncatt_get(fai_max,0,"institution")
datasource <- ncatt_get(fai_max,0,"source")
references <- ncatt_get(fai_max,0,"references")
history <- ncatt_get(fai_max,0,"history")
Conventions <- ncatt_get(fai_max,0,"Conventions")


library(lattice)
library(RColorBrewer)
library(CFtime)
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
cf
timestamps <- as_timestamp(cf) # get character-string times
class(timestamps)

time_cf <- parse_timestamps(cf, timestamps) # parse the string into date components
time_cf
class(time_cf)

fai_array[fai_array==fillvalue$value] <- NA
length(na.omit(as.vector(fai_array[,,1])))
