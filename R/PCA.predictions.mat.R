PCA.predictions.mat <-
function (X, scaled.mat = FALSE, e.vects = 1:2) 
{
    X <- as.matrix(X)
    unscaled.X <- X
    means <- apply(X, 2, mean)
    sd <- sqrt(apply(X, 2, var))
    if (scaled.mat) 
        X <- scale(X)
    else {
        X <- scale(X, scale = FALSE)
        sd <- rep(1, ncol(X))
    }
    X.svd <- svd(X)
    X.hat <- X %*% X.svd$v[, e.vects] %*% t(X.svd$v[, e.vects])
    sweep(sweep(X.hat, 2, sd, "*"), 2, means, "+")
}
