# 05/05/2020
# test ‘dirichletprocess’ package

X = c(runif(50, min = -0.5, max = 0), runif(50, min = 0.5, max = 1))
Y = c(runif(50, min = -0.5, max = 0), runif(50, min = 0.5, max = 1))
dat = data.frame(X=X, Y=Y)
datTrans = scale(dat)

dpCl = DirichletProcessMvnormal(datTrans)

dpCl = Fit(dpCl, 2000, progressBar = TRUE)
plot(dpCl)


faithfulTrans <- scale(faithful)
dpCluster <-  DirichletProcessMvnormal(faithfulTrans)
dpCluster <- Fit(dpCluster, 2000, progressBar = TRUE)
plot(dpCluster)
