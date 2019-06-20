SMACOF <-
function (DeltaMat, Wmat, X.ini, tol = 1e-07, iter = 10, scaled.mat = TRUE) 
{
    eta.sq.delta <- sum(DeltaMat * DeltaMat * Wmat)
    n <- nrow(DeltaMat)
    X <- X.ini
    if (scaled.mat) 
        X <- scale(X)
    V <- -(Wmat + t(Wmat)) + diag(apply((Wmat + t(Wmat)), 1, 
        sum))
    Vplus <- solve(V + matrix(1, nrow = n, ncol = n))
    W.delta.mat <- DeltaMat * Wmat
    tel <- 0
    stress.begin <- Inf
    repeat {
        tel <- tel + 1
        if (tel > iter) {
            cat("Maksimum aantal toegelate iterasies is bereik! \n")
            break
        }
        D.z <- sqrt(Pythagoras.dist(X))
        B.z.mat <- ifelse(upper.tri(W.delta.mat), W.delta.mat/D.z, 
            W.delta.mat)
        B.z.mat <- B.z.mat + t(B.z.mat)
        B.z.mat <- -B.z.mat + diag(apply(B.z.mat, 2, sum))
        eta.sq.x <- sum(diag(t(X) %*% V %*% X))
        rho.x <- sum(diag(t(X) %*% B.z.mat %*% X))
        sigma.r <- eta.sq.delta + eta.sq.x - 2 * rho.x
        X.up <- Vplus %*% B.z.mat %*% X
        if (scaled.mat) 
            X.up <- scale(X.up)
        X <- X.up
        if ((stress.begin - sigma.r) < tol) 
            break
        stress.begin <- sigma.r
    }
    list(B.z.mat = B.z.mat, X.up = X.up, sqrt(Pythagoras.dist(X.up)), 
        sigma.r = sigma.r, iterasies = tel)
}
