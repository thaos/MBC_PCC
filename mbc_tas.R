library(ncdf4)
library(rvinecopulib)
source("rvinefun.R")
chosenmonth = "12"

get_month <- function(time){
  strsplit(as.character(time), " ")  |> 
    lapply(FUN = function(x) x[[1]])  |> 
    unlist() |>
    strsplit(split = "-") |>
    lapply(FUN = function(x) x[[2]]) |>
    unlist() 
}

ncfiles <- c(
  "tas_day_CNRM-CM5_piControl_r1i1p1_19001999_europe.nc",
  "tas_day_IPSL-CM5A-MR_piControl_r1i1p1_18000101-18191231_europe.nc",
  "tas_day_IPSL-CM5A-MR_piControl_r1i1p1_20000101-20191231_coarsegrain_europe.nc",
  "tas_day_IPSL-CM5A-MR_piControl_r1i1p1_20000101-20191231_europe.nc"
)


nc <- nc_open(ncfiles[1])
nlon <- nlat <- 28
lon <- (ncvar_get(nc, "lon"))[1:nlon]
lat <- (ncvar_get(nc, "lat"))[1:nlat]
nc_close(nc)

ltas <- lapply(
  ncfiles,
  function(file, nlon, nlat){
    nc <- nc_open(file)
    month <- ncdf4.helpers::nc.get.time.series(nc) |>
      get_month()
    tas <- ncvar_get(nc, "tas")  |>
      aperm(perm = c(3, 2, 1)) |> 
      (function(arr){ matrix(arr[, 1:nlon, 1:nlat], nrow = nrow(arr)) })() |>
      (function(x) x[month == chosenmonth, ])()
    nc_close(nc)
    return(tas)
  },
  nlon = nlon, nlat = nlat
)

x1 <- ltas[[1]]
y1 <- ltas[[2]]
x2 <- ltas[[3]]
y2 <- ltas[[4]]

image.plot <- fields::image.plot
zlim <- range(x1, x2, y1, y2)
par(mfrow = c(1, 2))
image.plot(lon, lat, matrix(x1[100,], nlon, nlat), zlim = zlim)
lines(maps::map(plot = FALSE))
image.plot(lon, lat, matrix(y1[100,], nlon, nlat), zlim = zlim)
lines(maps::map(plot = FALSE))



ux1 <- apply(x1, 2, function(x) rank(x)/(length(x) + 1))
ux2 <- apply(x2, 2, function(x) rank(x)/(length(x) + 1))
uy1 <- apply(y1, 2, function(x) rank(x)/(length(x) + 1))
uy2 <- apply(y2, 2, function(x) rank(x)/(length(x) + 1))

rm(x1, x2, y1, y2)
gc(TRUE)

zlim <- c(0, 1)
par(mfrow = c(1, 2))
image.plot(lon, lat, matrix(ux1[1,], length(lon), length(lat)), zlim = zlim)
lines(maps::map(plot = FALSE))
image.plot(lon, lat, matrix(uy1[1,], length(lon), length(lat)), zlim = zlim)
lines(maps::map(plot = FALSE))

family_set <- c("indep", "gaussian", "t","clayton", "frank", "joe",
                "bb1", "bb6", "bb7", "bb8")
family_set <- "nonparametric"
par_method <- "mle"
d <- ncol(uy1)
vy1 <- vinecop(
  uy1[, 1:d],
  var_types = rep("c", d),
  family_set = family_set,
  structure = NA, # dvine_structure(1:d, trunc_lvl = Inf),
  par_method = par_method,
  trunc_lvl = NA,
  show_trace = TRUE,
)
# summary(vy1)
# print(vy1$structure)
# print(vy1$structure$order)
# get_m(vy1)

uy1s <- rvinecop(nrow(uy1), vy1)
pairs(uy1s[, 1:10])
# image.plot(lon, lat, matrix(uy1s[100, ], length(lon), length(lat)), zlim = zlim)
# lines(maps::map(plot = FALSE))


wy1s <- rosenblatt(uy1s[, 1:d], vy1)
pairs(wy1s[, 1:10])

wy1 <- rosenblatt(uy1[, 1:d], vy1)
pairs(wy1[1:100, 1:10])

uy1r <- inverse_rosenblatt(wy1s, vy1)
pairs(uy1r[, 1:10])

plot(uy1s[, vy1$structure$order[1]] - uy1r[,  vy1$structure$order[1]])
sd(uy1s[, vy1$structure$order[1]] - uy1r[,  vy1$structure$order[1]])


vx1 <- vinecop(
  ux1[, 1:d],
  var_types = rep("c", d),
  family_set = family_set,
  structure = vy1$structure,
  par_method = par_method,
  trunc_lvl = vy1$structure$trunc_lvl,
  show_trace = TRUE,
)

ux1s <- rvinecop(300, vx1)
wx1s <-rosenblatt(ux1s[, 1:d], vx1)
pairs(wx1s[1:300, 1:10])
ux1r <-inverse_rosenblatt(wx1s, vx1)
pairs(ux1r[, 1:10])
plot(ux1s[, vy1$structure$order[1]] - ux1r[,  vy1$structure$order[1]])

wx1 <- rosenblatt(ux1[, 1:d], vx1)
pairs(wx1[101:600, 1:10])
uy1r <- inverse_rosenblatt(wx1, vy1)
pairs(uy1r[, 1:10])
pairs(ux1[, 1:10])

zlim <- c(0, 1)
itime <- 1
par(mfrow = c(1, 3))
image.plot(lon, lat, matrix(ux1[itime, ], length(lon), length(lat)), zlim = zlim)
lines(maps::map(plot = FALSE))
image.plot(lon, lat, matrix(uy1[itime, ], length(lon), length(lat)), zlim = zlim)
lines(maps::map(plot = FALSE))
image.plot(lon, lat, matrix(uy1r[itime, ], length(lon), length(lat)), zlim = zlim)
lines(maps::map(plot = FALSE))

# print(mean((uy1r - uy1)^2))
# print(mean((ux1 - uy1)^2))

print(mean((cor(ux1) - cor(uy1))^2))
print(mean((cor(uy1r) - cor(uy1))^2))
print(mean((cor(uy1s) - cor(uy1))^2))

