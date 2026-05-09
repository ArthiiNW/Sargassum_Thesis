library(ncdf4)

fai<-nc_open("data/min_FAI_2020-26.nc")
print(fai)


lon <- ncvar_get(fai,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(fai,"latitude")
nlat <- dim(lat)
head(lat)

time <- ncvar_get(fai,"time")
tunits <- ncatt_get(fai,"time","units")
nt <- dim(time)

fai_array<-ncvar_get(fai,"nfai_min")
longname<- ncatt_get(fai, "nfai_min", "long_name")
dunits<-ncatt_get(fai, "nfai_min", "units")
fillvalue<-ncatt_get(fai, "nfai_min", "FillValue")

title <- ncatt_get(fai,0,"title")
institution <- ncatt_get(fai,0,"institution")
datasource <- ncatt_get(fai,0,"source")
references <- ncatt_get(fai,0,"references")
history <- ncatt_get(fai,0,"history")
Conventions <- ncatt_get(fai,0,"Conventions")


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

n_time <- fai$dim$time$len

max_fai     <- numeric(n_time)
n_positive  <- numeric(n_time)

for (t in seq_len(n_time)) {
  slice <- ncvar_get(fai, "nfai_min",
                     start = c(1, 1, t),
                     count = c(-1, -1, 1))
  max_fai[t]    <- suppressWarnings(max(slice, na.rm = TRUE))
  n_positive[t] <- sum(slice > 0, na.rm = TRUE)
}
nc_close(fai)

image(lon,lat,fai_slice, col=rev(brewer.pal(10,"RdBu")))
