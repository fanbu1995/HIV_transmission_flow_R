############################## -- 
## Generalized Spatial Fusion Model Framework for 
## Joint Analysis of Point and Areal Data, Appendix B

## Please cite: Wang, C., Puhan, M.A., and Furrer. R, (2017) Generalized Spatial Fusion Model Framework for Joint Analysis of Point and Areal Data, Spatial Statistics.

## Contains R code used to generate the simulated data

## Authors: Craig Wang, Milo Puhan, Reinhard Furrer (reinhard.furrer@math.uzh.ch)
############################## --
# n.area = number of areal regions
# dimension = dimension of the simulated region
# psill = partial sill of spatial process
# eff.range = effective range based on a exponential model
# nugget = 0, the nugget is incorporated into the response variance
# tau.sq = nugget
# beta = p by 1 matrix of coefficents for point response
# eta = p by 1 matrix of coefficients for logitstic regression of risk
# alpha = p by 1 matrix of coefficients for areal response, alpha2 is the coefficient of aggregated risk on areal response
# n.sample = number of locations
# n.neighbor = number of nearest neighbors
# L, L2 = the number of sampling points for the dataset dat.fusion and dat.fusionL2 (generate 2 datasets in one go for convenience)
library(spam)
library(deldir)
library(sp)
library(gstat)
library(ggplot2)
library(splancs)
library(boot)
library(maptools)
library(compare)
library(broom)
library(spdep)
library(spBayes)
library(fields)
library(rstan)
library(coda)
library(MASS)
library(Matrix)
# Main function -----------------------------------------------------------
# Simulate data for both fusion and non-fusion model according to parameters
simulate_data <- function(seed=2017, n.area=50, dimension=4000, psill=1, phi=900, nugget=0, tau.sq=1, L1=1, L2=5, tau_a.sq=0, sampling1="random", sampling2="random",
                         beta=rbind(1,5), alpha=rbind(1,0.5), alpha2=2, n.sample=500, n.pred=300, n.neighbor=5, eta=rbind(-2,5), eta2=-1, linkYQ=FALSE,f.cov="Exp",weight){
  
  # Supporting functions 
  voronoipolygons = function(centroids) {
    # generate spatialpolygon
    # function adopted from http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
    z = deldir(centroids[,1], centroids[,2])
    w = tile.list(z)
    
    for (i in 1:length(w)){
      w[[i]]$x[ w[[i]]$x < 0] <- 0 
      w[[i]]$y[ w[[i]]$y < 0] <- 0 
      w[[i]]$x[ w[[i]]$x > dimension] <- dimension 
      w[[i]]$y[ w[[i]]$y > dimension] <- dimension 
    }
    
    polys = vector(mode='list', length=length(w))
    for (i in seq(along=polys)) {
      pcrds = cbind(w[[i]]$x, w[[i]]$y)
      pcrds = rbind(pcrds, pcrds[1,])
      polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
    }
    SP = SpatialPolygons(polys)
    voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=centroids[,1], y=centroids[,2], row.names=sapply(slot(SP, 'polygons'), function(x) slot(x, 'ID'))))
    return(voronoi)
  }
  # get the neighborhood index for each location
  # s: location matrix, each row records coordinates of one point
  # m: number of neighborhood
  i_index <- function(i, s, n.neighbor) {
    #if (is.element("fields", installed.packages()[, 1]))
    #  library(fields)
    
    if(n.neighbor >= (i - 1)) {im <- 1:(i - 1)}
    else 	{
      dist <- rdist(s[c(1,i),], s[c(1:(i-1)), ])[-1,]
      im <- sort(order(dist)[1:n.neighbor])
    }
    return(im)
  }
  # distance matrix for location i and its neighbors 
  i_dist <- function(i, neighbor_index, s)	dist(s[c(i, neighbor_index[[i - 1]]), ])
  get_index_dist <- function(s, n.neighbor) {
    n = nrow(s)
    n.neighbor = min(n.neighbor, n-1)
    # get index of neighborhood
    neighbor_index <- sapply(2:n, i_index, s, n.neighbor)
    # get distance matrix for each i
    neighbor_dist <- sapply(2:n, i_dist, neighbor_index, s)
    return(list(i = neighbor_index, d = neighbor_dist))
  }
  # get the output for c function 
  get_neardistM <- function (ind, ind_distM_d) {
    if (ind < n.neighbor ){l = ind } else {l = n.neighbor}; 
    M_i <- rep(0, n.neighbor * (n.neighbor - 1) / 2);
    if (l == 1) {}
    else{
      M_i[1: (l * (l - 1) / 2)] <- 
        c(ind_distM_d[[ind]])[(l + 1): (l * (l + 1) / 2)]
    }
    return(M_i)
  }
  get_neardist <- function (ind, ind_distM_d) {
    if (ind < n.neighbor ){l = ind } else {l = n.neighbor}; 
    D_i <- rep(0, n.neighbor);
    D_i[1:l]<-c(ind_distM_d[[ind]])[1:l]
    return( D_i)
  }
  get_nearind <- function (ind, ind_distM_i) {
    if (ind < n.neighbor ){l = ind } else {l = n.neighbor}; 
    D_i <- rep(0, n.neighbor);
    D_i[1:l]<-c(ind_distM_i[[ind]])[1:l]
    return( D_i)
  }
  
  set.seed(seed)

  # generate latent spatial process 
  # resolution <- dimension*2/100 # result approximate 10,000 locations
  if (f.cov=="Mat"){
    vario.model <- vgm(psill=psill,model="Mat",range=phi/sqrt(3),nugget = 0, kappa = 1.5) # see plots.R in poster
  } else {
    vario.model <- vgm(psill=psill,model="Exp",range=phi,nugget = 0)
  }
#   xy <- data.frame(x=runif(min(n.sample*100, 100000+n.sample),0,dimension),
#                    y=runif(min(n.sample*100, 100000+n.sample),0,dimension))
  xy <- data.frame(x=runif(100000,0,dimension),
                   y=runif(100000,0,dimension))
  
  n.pts <- floor(sqrt(n.pred))
  xy.pred <- expand.grid(seq(dimension/(n.pts-1)/2, dimension-dimension/(n.pts-1)/2, length = n.pts),
                         seq(dimension/(n.pts-1)/2, dimension-dimension/(n.pts-1)/2, length = n.pts))
  names(xy.pred) <- c("x","y")
  xy <- rbind(xy, xy.pred)
  
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vario.model, nmax=50) # set to 50 to reduce computation time
  mrf <- predict(g.dummy, newdata=xy, nsim=1)
  
  # generate aggregated latent process -
  centroids.x <- runif(n.area,0,dimension)
  centroids.y <- runif(n.area,0,dimension)
  poly <- voronoipolygons(cbind(centroids.x,centroids.y))
  
  pips <- lapply(seq(1:n.area),function(x)  pip(xy[,c("x","y")],cbind(x=poly@polygons[[x]]@Polygons[[1]]@coords[,1],y=poly@polygons[[x]]@Polygons[[1]]@coords[,2])))
  pips.index <- lapply(seq(1:n.area),function(x) inpip(xy[,c("x","y")],cbind(x=poly@polygons[[x]]@Polygons[[1]]@coords[,1],y=poly@polygons[[x]]@Polygons[[1]]@coords[,2]),bound=T))
  
  aD_matrix <- matrix(0,n.area,nrow(xy))
  for (i in 1:n.area){
    aD_matrix[i,pips.index[[i]]] <- 1 
  }
  aD_matrix.full <-  diag(1/rowSums(aD_matrix)) %*% aD_matrix  # aggregation matrix that computes average mrf per area 
  mean.mrf <- aD_matrix.full %*% mrf$sim1
  
  poly.tidy <- tidy(poly)
  poly.tidy$wbar <- mean.mrf[as.numeric(poly.tidy$id)]
  
  # L1 sampling locations with a centroid plus (L1-1) points within each area
  if (L1==1){
    sample.locsL1 <- cbind(centroids.x,centroids.y)
    colnames(sample.locsL1) <- c("x","y")
  } else {
    sample.locsL1 <- data.frame(do.call(rbind,sapply(1:n.area, function(i) {
      polysample <- spsample(poly@polygons[[i]],n=L1-1,type = sampling1)@coords
      cbind(id=i,rbind(polysample,c(centroids.x[i],centroids.y[i])),realL=nrow(polysample)+1)}
      ,simplify = F)))
    colnames(sample.locsL1)[2:3] <- c("x","y")
  }
  
  sample.locsL2 <- data.frame(do.call(rbind,sapply(1:n.area, function(i) {
    polysample <- spsample(poly@polygons[[i]],n=L2-1,type = sampling2)@coords
    cbind(id=i,rbind(polysample,c(centroids.x[i],centroids.y[i])),realL=nrow(polysample)+1)}
    ,simplify = F)))
  colnames(sample.locsL2)[2:3] <- c("x","y")
  sample.ind <- sample(nrow(xy),n.sample)
  
  # generate areal response 
  Xn <- cbind(rep(1,nrow(xy)),rnorm(nrow(xy),0,1))
  w <- mrf$sim1
  y <- rnorm(nrow(xy),Xn %*% beta + w, sqrt(tau.sq))
  ind.risk <- inv.logit(Xn %*% eta + eta2 * y)
  D <- rbinom(nrow(xy),1,ind.risk)
  
  
  Q <- NULL
  Xm <- cbind(rep(1,n.area),rnorm(n.area,0,1))
  
  for (i in 1:n.area){
    if (linkYQ){
      Q[i]<-rpois(1,exp((Xm %*% alpha)[i] + (alpha2 * aD_matrix.full %*% ind.risk)[i] + poly.tidy$wbar[as.numeric(poly.tidy$id)==i][1]))
    }
    else {Q[i]<-rpois(1,exp((Xm %*% alpha)[i] + poly.tidy$wbar[as.numeric(poly.tidy$id)==i][1] + rnorm(1,0,sd=sqrt(tau_a.sq)) ))}
  }
  poly.tidy$Q <- Q[as.numeric(poly.tidy$id)]
  
  
  if (L1==1){
    A0L1 <- cbind(diag(n.sample),matrix(0,n.sample,n.area*L1))
    A1L1.2nd <- as.matrix(bdiag(lapply(1:n.area, function(i) t(matrix(rep(1/L1,L1)))))) 
    A1L1 <- cbind(matrix(0,n.area,n.sample),A1L1.2nd)
  } else {
    idrealL1 <- aggregate(realL~id,data=sample.locsL1,length) # number of location per polygon
    A0L1 <- cbind(diag(n.sample),matrix(0,n.sample,nrow(sample.locsL1)))
    A1L1.2nd <- as.matrix(bdiag(lapply(1:n.area, function(i) t(matrix(rep(1/idrealL1$realL[i],idrealL1$realL[i]))))))  # the non-zero part of A1, block diag
    A1L1 <- cbind(matrix(0,n.area,n.sample),A1L1.2nd)
  }
  
  idrealL2 <- aggregate(realL~id,data=sample.locsL2,length)
  A0L2 <- cbind(diag(n.sample),matrix(0,n.sample,nrow(sample.locsL2))) # number of location per polygon
  A1L2.2nd <- as.matrix(bdiag(lapply(1:n.area, function(i) t(matrix(rep(1/idrealL2$realL[i],idrealL2$realL[i]))))))  # the non-zero part of A1, block diag
  A1L2 <- cbind(matrix(0,n.area,n.sample),A1L2.2nd)
  
  fit <- lm(y[sample.ind]~Xn[sample.ind,]-1)
  beta.starting <- coefficients(fit)
  alpha.starting <- coefficients(glm(Q~Xm-1,family="poisson"))
  beta.tuning <- t(chol(vcov(fit)))
  
  # NNGP 
  ind_distM <- get_index_dist(mrf[sample.ind,c("x","y")], n.neighbor)
  neardistM <- sapply(1: (n.sample- 1), get_neardistM, ind_distM$d)
  neardist <- sapply(1:(n.sample-1), get_neardist, ind_distM$d)
  nearind <- sapply(1: (n.sample-1), get_nearind, ind_distM$i)
  
  FCOV = ifelse(f.cov=="Exp",1,2)
  
  dat <- list(N = n.sample, M = n.neighbor, P = 2, Y = y[sample.ind], X = Xn[sample.ind,1:2],
              nearind = t(nearind), neardist = t(neardist), neardistM = t(neardistM), FCOV = FCOV)
  
  ind_distML1 <- get_index_dist(rbind(mrf[sample.ind,c("x","y")],sample.locsL1[,c("x","y")]), n.neighbor)
  neardistML1 <- sapply(1: (n.sample++nrow(sample.locsL1) - 1), get_neardistM, ind_distML1$d)
  neardistL1 <- sapply(1:(n.sample++nrow(sample.locsL1)-1), get_neardist, ind_distML1$d)
  nearindL1 <- sapply(1: (n.sample++nrow(sample.locsL1)-1), get_nearind, ind_distML1$i)
  
  ind_distML2 <- get_index_dist(rbind(mrf[sample.ind,c("x","y")],sample.locsL2[,c("x","y")]), n.neighbor)
  neardistML2 <- sapply(1: (n.sample+nrow(sample.locsL2) - 1), get_neardistM, ind_distML2$d)
  neardistL2 <- sapply(1:(n.sample+nrow(sample.locsL2)-1), get_neardist, ind_distML2$d)
  nearindL2 <- sapply(1: (n.sample+nrow(sample.locsL2)-1), get_nearind, ind_distML2$i)
  
  ind_distM_sample <- get_index_dist(sample.locsL2[,c("x","y")], n.neighbor)
  neardistM_sample <- sapply(1: (nrow(sample.locsL2) - 1), get_neardistM, ind_distM_sample$d)
  neardist_sample <- sapply(1:(nrow(sample.locsL2)-1), get_neardist, ind_distM_sample$d)
  nearind_sample <- sapply(1: (nrow(sample.locsL2)-1), get_nearind, ind_distM_sample$i)
  
  aD_matrix.partial <- diag(1/rowSums(aD_matrix[,sample.ind])) %*% aD_matrix[,sample.ind] # aggregation matrix based on the observations
  aD_matrix.partial[is.nan(aD_matrix.partial)]<-0
  
  temp <- over(SpatialPoints(mrf[sample.ind,c("x","y")]), SpatialPolygons(poly@polygons))
  weights <- 1/table(temp)
  weight <- sapply(temp, function(x) weights[as.numeric(names(weights))==x])
  weight <- c(weight, rep(1/L2, n.area*L2))
  if (linkYQ){
    dat.fusion <- list(n = n.sample, M = n.neighbor, a = n.area, aL = n.area*L, pn = 2, pa = 2, Y = y[sample.ind], Xn = Xn[sample.ind,1:2], Xa = Xm, D = D[sample.ind],
                       nearind = t(nearind), neardist = t(neardist), neardistM = t(neardistM), Q = Q, A0 = A0, A1=A1, aD_matrix = aD_matrix.partial)
  } else {  dat.fusionL1 <- list(n = n.sample, M = n.neighbor, a = n.area, aL = nrow(sample.locsL1), pn = 2, pa = 2, Y = y[sample.ind], Xn = Xn[sample.ind,1:2], Xa = Xm,
                               nearind = t(nearindL1), neardist = t(neardistL1), neardistM = t(neardistML1), Q = Q, A0 = A0L1, A1=A1L1, aD_matrix = aD_matrix.partial, FCOV=FCOV,
                               samplelocs = sample.locsL1)
  }
  dat.fusionL2 <- list(n = n.sample, M = n.neighbor, a = n.area, aL = nrow(sample.locsL2), pn = 2, pa = 2, Y = y[sample.ind], Xn = Xn[sample.ind,1:2], Xa = Xm,
                       nearind = t(nearindL2), neardist = t(neardistL2), neardistM = t(neardistML2), 
                       nearind_sample = t(nearind_sample), neardist_sample = t(neardist_sample), neardistM_sample = t(neardistM_sample),
                       Q = Q, A0 = A0L2, A1=A1L2, aD_matrix = aD_matrix.partial, true_tau2 = tau.sq, true_sigma2= psill, FCOV=FCOV, samplelocs = sample.locsL2, weight=weight)
  return(list(dat=dat, dat.fusionL1 = dat.fusionL1, dat.fusionL2 = dat.fusionL2, sample.ind = sample.ind, mrf = mrf, y = y, X = Xn,
              beta.starting = beta.starting, alpha.starting = alpha.starting, beta.tuning = beta.tuning, poly = poly, mean.w = mean.mrf, 
              pred.loc = xy.pred, pred.ind = seq(nrow(xy) - nrow(xy.pred) + 1,  nrow(xy))))
}

# iter.stan <- 2000
# simulation <- simulate_data(n.sample=100)
# create.starting <- function(alpha.starting, beta.starting, tau.sq, psill, phi){
#   if (tau.sq == 0){tau.sq <- 0.01}
#   alpha.starting <- abs(alpha.starting)
#   beta.starting <- abs(beta.starting)
#   starting <- replicate(4,list("alpha"=rnorm(2,alpha.starting,alpha.starting/20),"beta"=rnorm(2,beta.starting,beta.starting/20), 
#                                               "tau_sq" = (rnorm(1,tau.sq,tau.sq/20)), "sigma_sq"=(rnorm(1,psill,psill/20)), 
#                                               "phi"=rnorm(1,phi,phi/20)), simplify = FALSE)
#   return(starting)
# }
# model.fusion <- stan_model(file = "fusion.stan", model_name = "fusion", verbose = T)
# starting.fusion <- create.starting(simulation$alpha.starting, simulation$beta.starting, 1, 1, 300)
# fit.fusion <- sampling(model.fusion, pars = c("beta","alpha","phi","tau_sq","sigma_sq","w"),
#                        data = simulation$dat.fusionL2, init = starting.fusion, seed=2016, chains = 4, cores = 4, 
#                        iter = iter.stan, control = list(stepsize = 0.01, max_treedepth = 12, adapt_delta = 0.95), verbose = TRUE) 
# monitor(extract(fit.fusion, pars=c("alpha","beta","phi","tau_sq","sigma_sq"), permute=FALSE))
# prediction <- short_pred.fusion.w(fit.fusion, simulation)
# plot(simulation$mrf$sim1[simulation$pred.ind], apply(prediction$w, 1,median))
