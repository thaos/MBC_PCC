library(rvinecopulib)
source("rvinefun.R")
## simulate dummy data
d <- 10
nsim <- 100
w <- matrix(runif(nsim * d), nsim, d)

x <- rnorm(30) * matrix(1, 30, d) + 0.5 * matrix(rnorm(30 * d), 30, d)
u <- pseudo_obs(x)
pairs(u)

## fit a C-vine model 
vc <- vinecop(u, family = "clayton", structure = cvine_structure(1:d, trunc_lvl = Inf)
)
summary(vc)
mstruct <- as_rvine_matrix(vc$structure)
mstruct
nsim <- 1000
u1 <- rvinecop(nsim, vc)
pairs(u1)


u2 <- w2u_cvine(w, vc)
pairs(u2)
vc1 <- vinecop(u1, family = "clayton", structure = cvine_structure(1:d, trunc_lvl = Inf)
)
vc2 <- vinecop(u2, family = "clayton", structure = cvine_structure(1:d, trunc_lvl = Inf)
)
summary(vc)
summary(vc1)
summary(vc2)



w2 <- u2w_cvine(u2, vc)
pairs(w2)
plot(w2[, vc$structure$order[1]] - w[,  vc$structure$order[1]])

bicop <- vc$pair_copulas[[d - 1]][[1]]
fam <- bicop$family
params <- bicop$parameters
rot <- bicop$rotation
test = hbicop(
  cbind(
    hbicop(
      cbind(
        u2[, 2],
        u2[, 3]
      ), 
      2, fam, rot, params
    ),
    u2[,3]
  ),
  2, fam, rot, params, inverse = TRUE
)


## fit a D-vine model 
vc <- vinecop(u, family = "clayton", structure = dvine_structure(1:d, trunc_lvl = Inf)
)
summary(vc)
mstruct <- as_rvine_matrix(vc$structure)
mstruct
u1 <- rvinecop(nsim, vc)
pairs(u1)



u2 <- w2u_dvine(w, vc)
pairs(u2)

vc1 <- vinecop(u1, family = "clayton", structure = dvine_structure(1:d, trunc_lvl = Inf)
)
vc2 <- vinecop(u2, family = "clayton", structure = dvine_structure(1:d, trunc_lvl = Inf)
)
summary(vc)
summary(vc1)
summary(vc2)



w2 <- u2w_dvine(u2, vc)
pairs(w2)
plot(w2[, vc$structure$order[1]] - w[,  vc$structure$order[1]])






## fit a R-vine model 
vc <- vinecop(u, family = "clayton")
summary(vc)
mstruct <- as_rvine_matrix(vc$structure)
mstruct
nsim <- 100
u1 <- rvinecop(nsim, vc)
pairs(u1)




u2 <- w2u_rvine(w, vc)
pairs(u2)

vc1 <- vinecop(u1, family = "clayton", structure = vc$structure)
vc2 <- vinecop(u2, family = "clayton", structure = vc$structure)
summary(vc)
summary(vc1)
summary(vc2)


w2 <- u2w_rvine(u2, vc)
pairs(w2)
plot(w2[, vc$structure$order[1]] - w[,  vc$structure$order[1]])




