BetweenWithinSS <-
function (X, G) 
{
    X.center <- scale(X, center = TRUE, scale = FALSE)
    n <- nrow(X)
    g <- ncol(G)
    ng <- apply(G, 2, sum)
    SS11 <- t(G) %*% G
    SS12 <- t(G) %*% X.center
    SS.T <- t(X) %*% X.center
    SS.B <- t(SS12) %*% solve(SS11) %*% SS12
    SS.W <- SS.T - SS.B
    MeansMat.Xcenter <- solve(diag(ng)) %*% t(G) %*% X.center
    MeansMat <- solve(diag(ng)) %*% t(G) %*% as.matrix(X)
    list(SS.T = SS.T, SS.B = SS.B, SS.W = SS.W, MeansMat = MeansMat, 
        MeansMat.Xcenter = MeansMat.Xcenter)
}
