
w2u_cvine <- function(w, cvine_fit) {
  d <- ncol(w)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- cvine_fit$structure$order
  revordercol <- order(ordercol)
  v <- array(NA, dim = c(d, d, nrow(w)))
  v[1,  1,] <-  w[, ordercol][, revcol[1]]
  for( i in 2:d ){
    v[i, i, ] <-  w[, ordercol][, revcol[i]]
    for(k in (i-1):1){
      bicop <- vc$pair_copulas[[k]][[revcol[i]]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      v[k, i, ] <- hbicop(
        cbind(
          v[k + 1, i, ],
          v[k, k, ]
        ),
        2, fam, rot, params, inverse = TRUE
      )
    }
  }
  u <- t(v[1,,])[, revcol][, revordercol]
  return(u)
}

u2w_cvine <- function(u, cvine_fit) {
  d <- ncol(u)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- cvine_fit$structure$order
  revordercol <- order(ordercol)
  v <- array(NA, dim = c(d, d, nrow(w)))
  v[1,  ,] <-  t(u[, ordercol][, revcol])
  for( i in 2:d ){
    for ( k in  1:(i-1) ){
      bicop <- cvine_fit$pair_copulas[[k]][[ revcol[i] ]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      v[k + 1, i, ] <- hbicop(
        cbind(
          v[k, i, ],
          v[k, k, ]
        ),
        2, fam, rot, params
      )
    }
  }
  w <- mapply(function(i, j) v[i, j, ], 1:d, 1:d)[, revcol][, revordercol]
  return(w)
}

w2u_dvine <- function(w, dvine_fit) {
  d <- ncol(w)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- dvine_fit$structure$order
  revordercol <- order(ordercol)
  
  v <- v2 <- array(NA, dim = c(d, d, nrow(w)))
  v[1,  1,] <-  v2[1,  1, ] <-  w[, ordercol][, revcol[1]]
  for( i in 2:d ){
    v[i, i, ] <-  w[, ordercol][, revcol[i]]
    for(k in (i-1):1){
      bicop <- dvine_fit$pair_copulas[[k]][[revcol[i]]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      v[k, i, ] <- hbicop(
        cbind(
          v[k + 1, i, ],
          v2[k, i - 1, ]
        ),
        2, fam, rot, params, inverse = TRUE
      )
      if (i < d) {
        v2[k + 1, i, ] <- hbicop(
          cbind(
            v[k, i, ],
            v2[k, i - 1, ]
          ),
          1, fam, rot, params
        )
      }
    }
    v2[1, i, ] <-  v[1, i, ]
  }
  u <- t(v[1,,])[, revcol] [, revordercol]
  return(u)
}


u2w_dvine <- function(u, dvine_fit) {
  d <- ncol(u)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- dvine_fit$structure$order
  revordercol <- order(ordercol)
  
  v <- v2 <- array(NA, dim = c(d, d, nrow(u)))
  v[1,  ,] <-  v2[1,  , ] <-  t(u[, ordercol][, revcol])
  
  for( i in 2:d ){
    
    for ( k in  1:(i-1) ){
      
      bicop <- dvine_fit$pair_copulas[[k]][[revcol[i]]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      
      v[k + 1, i, ] <- hbicop(
        cbind(
          v[k, i, ],
          v2[k, i - 1, ]
        ),
        2, fam, rot, params
      )
      
      v2[k + 1, i, ] <- hbicop(
        cbind(
          v[k, i, ],
          v2[k, i - 1, ]
        ),
        1, fam, rot, params
      )
      
    }
  }  
  w <- mapply(function(i, j) v[i, j, ], 1:d, 1:d)[, revcol][, revordercol]
  return(w)
}


get_m <- function(vc) {
  struct <- vc$structure
  d <- struct$d
  m <- matrix(NA, d, d)
  #convtab <- (1:d)[struct$order]
  convtab <- d:1
  for(i in 1:length(struct$struct_array)){
    m[i,(i + 1):d] <- rev(struct$struct_array[[i]])
  }
  diag(m) <- d:1
  m[] <- convtab[m]
  return(m)
}
get_mtilde <- function(m) apply(m, 2, cummax)

w2u_rvine <- function(w, rvine_fit) {
  m <- get_m(rvine_fit)
  mtilde <- get_mtilde(m)
  d <- ncol(w)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- rvine_fit$structure$order
  revordercol <- order(ordercol)
  trunc_lvl <- rvine_fit$structure$trunc_lvl
  v <- v2 <- array(NA, dim = c(min(d, trunc_lvl), d, nrow(w)))
  v[1,  1,] <-  v2[1,  1, ] <-  w[, ordercol][, revcol[1]]
  for ( i in 2:d ) {
    v[min(i, trunc_lvl), i, ] <-  w[, ordercol][,  revcol[i]]
    for ( k in min(i - 1, trunc_lvl - 1):1 ) {
      bicop <- rvine_fit$pair_copulas[[ k ]][[ revcol[i] ]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      if ( m[k, i] == mtilde[k, i] ) {
        v[k, i, ] <- hbicop(
          cbind(
            v[k + 1, i, ],
            v[k, mtilde[k, i], ]
          ),
          2, fam, rot, params, inverse = TRUE
        )
      } else {
        v[k, i, ] <- hbicop(
          cbind(
            v[k + 1, i, ],
            v2[k, mtilde[k, i], ]
          ),
          2, fam, rot, params, inverse = TRUE
        )
      }
      if ( any(is.nan(v[k, i,])) ){
        stop(paste0("NaN in v: i = ", i, ", k = ", k ))
      }
      if ( i < d ) {
        if (m[k, i] == mtilde[k, i]){
          v2[k + 1, i, ] <- hbicop(
            cbind(
              v[k, i, ],
              v[k, mtilde[k, i], ]
            ),
            1, fam, rot, params
          )
        } else {
          v2[k + 1, i, ] <- hbicop(
            cbind(
              v[k, i, ],
              v2[k, mtilde[k, i], ]
            ),
            1, fam, rot, params
          )
        }
        if ( any(is.nan(v2[k + 1, i,])) ) {
          stop(paste0("NaN in v2: i = ", i, ", k = ", k ))
        }
      }
      
    }
    v2[1, i, ] <-  v[1, i, ]
    
  }
  # reorder columns is in w
  u <- t(v[1, revcol, ]) [, revordercol]
  return(u)
}


u2w_rvine <- function(u, rvine_fit) {
  m <- get_m(rvine_fit)
  mtilde <- get_mtilde(m)
  d <- ncol(u)
  # in rvinecopulib, first variable has index 4 
  # while in Claudia Czado's book first variable has index 1
  # we follow the notation/indexes of  Claudia Czado's book
  revcol <- d:1
  # rvinecopulib permutes the columns
  ordercol <- rvine_fit$structure$order
  revordercol <- order(ordercol)
  trunc_lvl <- rvine_fit$structure$trunc_lvl
  v <- v2 <-  array(NA, dim = c(min(d, trunc_lvl), d, nrow(u)))
  v[1,  ,] <- v2[1,  ,] <- t(u[, ordercol][, revcol])
  for ( i in 2:d ) {
    for ( k in  1:min(i - 1, trunc_lvl - 1) ) {
      bicop <- rvine_fit$pair_copulas[[k]][[ revcol[i] ]]
      fam <- bicop$family
      params <- bicop$parameters
      rot <- bicop$rotation
      if ( m[k, i] == mtilde[k, i] ){
        v[k + 1, i, ] <- hbicop(
          cbind(
            v[k, i, ],
            v[k, mtilde[k, i], ]
          ),
          2, fam, rot, params
        )
      } else {
        v[k + 1, i, ] <- hbicop(
          cbind(
            v[k, i, ],
            v2[k, mtilde[k, i], ]
          ),
          2, fam, rot, params
        )
      }
      if (m[k, i] == mtilde[k, i]){
        v2[k + 1, i, ] <- hbicop(
          cbind(
            v[k, i, ],
            v[k, mtilde[k, i], ]         
          ),
          1, fam, rot, params
        )
      } else {
        v2[k + 1, i, ] <- hbicop(
          cbind(
            v[k, i, ],
            v2[k, mtilde[k, i], ]
          ),
          1, fam, rot, params
        )
      }
    }
  }
  w <- mapply(function(i, j) v[i, j, ], pmin(1:d, trunc_lvl), 1:d)[, revcol][, revordercol]
  return(w)
}