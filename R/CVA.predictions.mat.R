CVA.predictions.mat <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    weightedCVA = c("weighted", "unweightedI", "unweightedCent")) 
{
    weightedCVA <- weightedCVA[1]
    X <- as.matrix(X)
    unscaled.X <- X
    X.cent <- scale(X, center = TRUE, scale = FALSE)
    X.means <- apply(X, 2, mean)
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    K <- min(p, J - 1)
    Gmat <- G
    Nmat <- t(Gmat) %*% Gmat
    XcentBar.groups <- solve(Nmat) %*% t(Gmat) %*% X.cent
    SSP.T <- t(X.cent) %*% X.cent
    SSP.B <- t(XcentBar.groups) %*% Nmat %*% XcentBar.groups
    SSP.W <- SSP.T - SSP.B
    Wmat <- SSP.W
    svd.Wmat <- svd(Wmat)
    lambdamatI <- diag(svd.Wmat$d)
    Lmat <- svd.Wmat$u %*% solve(sqrt(lambdamatI))
    if (weightedCVA == "weighted") {
        Cmat <- Nmat
        Csqrt <- sqrt(Nmat)
    }
    if (weightedCVA == "unweightedI") {
        Cmat <- diag(J)
        Csqrt <- Cmat
    }
    if (weightedCVA == "unweightedCent") {
        Cmat <- diag(J) - matrix(1/J, nrow = J, ncol = J)
        Csqrt <- Cmat
    }
    if (is.na(match(weightedCVA, c("weighted", "unweightedI", 
        "unweightedCent")))) 
        stop(" Argument 'weightedCVA' must be one of 'weighted','unweightedI','unweightedCent' ")
    svd.step2 <- svd(t(Lmat) %*% t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups %*% Lmat)
    Vmat <- svd.step2$v
    lambdamat <- diag(svd.step2$d)
    svd.2sided <- Eigen.twosided(t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups, Wmat)
    Mmat <- svd.2sided$W
    lambdamat.2sided <- svd.2sided$Lambda.mat
    predictions.for.XcentBar.groups <- vector("list", p)
    predictions.for.centredX <- vector("list", p)
    predictions.for.originalX <- vector("list", p)
    reconstr.error.X <- rep(NA, p)
    reconstr.error.Xbar <- rep(NA, p)
    for (r in 1:p) {
        vec.temp <- rep(0, p)
        vec.temp[1:r] <- 1
        Jmat <- diag(vec.temp)
        XcentBarHat <- XcentBar.groups %*% Mmat %*% Jmat %*% 
            solve(Mmat)
        predictions.for.XcentBar.groups[[r]] <- XcentBarHat
        predict.centredX <- X.cent %*% Mmat %*% Jmat %*% solve(Mmat)
        predictions.for.centredX[[r]] <- predict.centredX
        predictions.for.originalX[[r]] <- scale(predict.centredX, 
            center = -attr(X.cent, "scaled:center"), scale = FALSE)
        dimnames(predictions.for.centredX[[r]])[[2]] <- dimnames(X)[[2]]
        dimnames(predictions.for.originalX[[r]])[[2]] <- dimnames(X)[[2]]
        reconstr.error.X[r] <- sum((X.cent - predict.centredX)^2)
        reconstr.error.Xbar[r] <- sum((XcentBar.groups - XcentBarHat)^2)
    }
    list(predictions.for.XcentBar.groups = predictions.for.XcentBar.groups, 
        predictions.for.centredX = predictions.for.centredX, 
        predictions.for.originalX = predictions.for.originalX, 
        reconstr.error.Xbar = reconstr.error.Xbar, reconstr.error.X = reconstr.error.X)
}
