############################## -- 
## Generalized Spatial Fusion Model Framework for 
## Joint Analysis of Point and Areal Data, Appendix B

## Please cite: Wang, C., Puhan, M.A., and Furrer. R, (2017) Generalized Spatial Fusion Model Framework for Joint Analysis of Point and Areal Data, Spatial Statistics.

## Contains R code used to generate posterior distributions of 
## the latent spatial process, and responses on predicted locations

## Authors: Craig Wang, Milo Puhan, Reinhard Furrer (reinhard.furrer@math.uzh.ch).	
############################## --
# tau.sq = nugget
# beta = p by 1 matrix of coefficents for point response
# n.sample = number of locations
# n.neighbor = number of nearest neighbors

# prediction on latent process --------------------------------------------
# sigma_sq,tau_sq,phi are posterior samples
# pred is locations for be predicted
# locs is the reference set
# n.neighbor is m used
# iter is number of posterior samples
# w is a nrow(locs) by iter matrix

# prediction for nngp
predict_nngp_w <- function(sigma_sq, phi, pred, locs, iter, n.neighbor, w){
  # Helper functions
  get_index_dist <- function(pred, locs , n.neighbor, n) {
    i_index <- function(i, pred, locs , n.neighbor) {
      #if (is.element("fields", installed.packages()[, 1]))
      #  library(fields)
      dist <- rdist(pred[i,], locs)
      im <- sort(order(dist)[1:n.neighbor])
      return(im)
    }
    # distance matrix for location i and its neighbors 
    i_dist <- function(i, neighbor_index, pred, locs )	rdist(pred[i,], locs[neighbor_index[,i],])
    # get index of neighborhood
    neighbor_index <- sapply(1:n, i_index, pred, locs , n.neighbor)
    # get distance matrix for each i
    neighbor_dist <- sapply(1:n, i_dist, neighbor_index, pred, locs )
    return(list(i = neighbor_index, d = neighbor_dist))
  }
  get_neardistM <- function (ind, neighbor_index, locs) {
    M_i <- c(dist(locs[neighbor_index[,ind],])) # unlist the distance matrix of neighbor - neighbors
    return(M_i)
  }
  
  n <- nrow(pred)
  w.pred <- matrix(0,n,iter) # n_pred by iteration matrix
  ind_distM <- get_index_dist(pred, locs, n.neighbor, n)
  neardistM <- t(sapply(1:n, get_neardistM, ind_distM$i, locs))
  neardist <- t(ind_distM$d) 
  neardist[neardist==0] <- min(neardist[neardist!=0])/2 # set coincided location to half of minimum distance
  nearind <- t(ind_distM$i)
  Ft <- BtW <- numeric(iter)
  for (i in 1:n) { # for each predicted location
      C_nei <- diag(1.001,n.neighbor) # m by m
      h = 0;
      for (j in 1:(n.neighbor-1)){
        for (k in (j+1):n.neighbor){
          h = h + 1
          C_nei[j, k] <- exp(- neardistM[i, h])
          C_nei[k, j] <- C_nei[j, k]
        }
      }
      C_site_nei <- exp(- neardist[i, ])
      
      for (it in 1:iter){ # for each iteration
        C_site_nei_C_nei_inv <- solve(C_nei^(1/phi[it]), C_site_nei^(1/phi[it]))# m by m times m by n
        Ft[it] <- C_site_nei_C_nei_inv %*% C_site_nei^(1/phi[it])
        BtW[it] <- C_site_nei_C_nei_inv %*% w[nearind[i,],it] # 1 by m, m by 1
        }
      w.pred[i,] <- rnorm(it,BtW,sqrt(sigma_sq * (1-Ft)))
  }
  return(w.pred)
}

# fitted response p(y|theta) p(theta|Data)
generate_y <- function(w, tau.sq, beta, X){
  # w[,i] and beta[i,] represents the ith iteration for n locations
  n <- nrow(X)
  iter <- ncol(w)
  y.pred <- matrix(0,n,iter) # n_pred by iteration matrix
  for (i in 1:iter){
    y.pred[,i] <- rnorm(n, X %*% t(beta[i,]) + w[,i], sqrt(tau.sq[i]))
  }
  return(y.pred)
}

short_pred.fusion.w <- function(model, simulation){
  # predict with both w at observated locations and sampling locations
  fit.fusion.extract <- data.frame(extract(model))[,c("sigma_sq","tau_sq","phi")]
  fit.fusion.beta <- data.frame(extract(model,c("beta")))
  fit.fusion.w <- t(as.matrix(data.frame(extract(model,"w")$w)))
  pred.ind <- simulation$pred.ind
  pred.fusion.w <- predict_nngp_w(fit.fusion.extract$sigma_sq, fit.fusion.extract$phi, simulation$mrf[pred.ind, c("x","y")], 
                                  rbind(simulation$mrf[simulation$sample.ind, c("x","y")],simulation$dat.fusionL2$samplelocs[,c("x","y")]),
                                        iter = iter.stan*2, n.neighbor = simulation$dat$M,
                                  fit.fusion.w)
  pred.fusion.y <- generate_y(pred.fusion.w, fit.fusion.extract$tau_sq, fit.fusion.beta, simulation$X[simulation$pred.ind,])
  
  return(list(w = pred.fusion.w, y = pred.fusion.y))
}
