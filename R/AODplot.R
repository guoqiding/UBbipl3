AODplot <-
function (X = scale(Pine.data[, -1]), group.vec = Pine.data[, 
    1], X.new.samples = scale(Pine.data[, -1]), dim.biplot = 2, 
    dist = c("Pythagoras", "Clark", "SqrtL1"), e.vects = 1:ncol(X), 
    exp.factor = 1.2, label = T, label.size = 0.8, pch.samples = rep(0, 
        length(unique(group.vec))), pch.samples.col = c("red", 
        "green", "blue", "yellow", "purple"), pch.samples.size = 1, 
    pch.means = rep(15, length(unique(group.vec))), pch.means.col = c("red", 
        "green", "blue", "yellow", "purple"), pch.means.size = 1.5, 
    pos.labs = 1, Title = NULL, weight = c("unweighted", "weighted")) 
{
    dist <- dist[1]
    weight <- weight[1]
    if (dist == "Clark") {
        if (any(X <= 0)) 
            stop("Clark distance only defined for positive X values")
        if (any(X.new.samples <= 0)) 
            stop("Clark distance only defined for positive X.new.samples values")
    }
    if (!is.null(X.new.samples)) {
        if (is.vector(X.new.samples)) 
            X.new.samples <- matrix(X.new.samples, nrow = 1)
        if (ncol(X.new.samples) != ncol(X)) 
            stop("X.new.samples must have same number of columns than X")
    }
    if (dist == "Clark") 
        dist.expr <- function(xik, xjk) {
            ((xik - xjk)/(xik + xjk))^2
        }
    if (dist == "Pythagoras") 
        dist.expr <- function(xik, xjk) {
            (xik - xjk)^2
        }
    if (dist == "SqrtL1") 
        dist.expr <- function(xik, xjk) {
            abs(xik - xjk)
        }
    D.fun <- function(mat) {
        n <- nrow(mat)
        D <- matrix(0, nrow = n, ncol = n)
        for (i in 1:(n - 1)) for (j in (i + 1):n) D[i, j] <- -0.5 * 
            sum(dist.expr(mat[i, ], mat[j, ]))
        D + t(D)
    }
    dvec.fun <- function(vec, mat) {
        if (length(vec) != ncol(mat)) 
            stop("Incorrect sized vector")
        n <- nrow(mat)
        outvec <- rep(NA, n)
        for (i in 1:n) outvec[i] <- -0.5 * sum(dist.expr(vec, 
            mat[i, ]))
        outvec
    }
    Gmat <- indmat(group.vec)
    groupnames <- colnames(Gmat)
    samples.colvec <- apply(Gmat, 1, function(x) (1:ncol(Gmat))[x == 
        1])
    K <- ncol(Gmat)
    Nmat <- t(Gmat) %*% Gmat
    nvec <- diag(Nmat)
    Dmat <- D.fun(X)
    Dbar <- solve(Nmat) %*% t(Gmat) %*% Dmat %*% Gmat %*% solve(Nmat)
    Deltabar <- Dbar - 0.5 * (diag(diag(Dbar)) %*% matrix(1, 
        nrow = K, ncol = K) + matrix(1, nrow = K, ncol = K) %*% 
        diag(diag(Dbar)))
    T.ss <- -sum(Dmat)/sum(nvec)
    B.ss <- -(t(nvec) %*% Deltabar %*% nvec)/sum(nvec)
    W.ss <- T.ss - B.ss
    DoubCentDeltabar <- (diag(K) - matrix(1, ncol = K, nrow = K)/K) %*% 
        Deltabar %*% (diag(K) - matrix(1, ncol = K, nrow = K)/K)
    temp <- matrix(1, ncol = 1, nrow = K) %*% matrix(nvec, ncol = K, 
        nrow = 1)
    DoubCentDeltabar.weighted <- (diag(K) - temp/sum(nvec)) %*% 
        Deltabar %*% (diag(K) - temp/sum(nvec))
    PCO.unweighted <- svd(DoubCentDeltabar)
    PCO.weighted <- svd(DoubCentDeltabar.weighted)
    Lambdamat.unweighted <- diag(PCO.unweighted$d[-K])
    Lambdamat.weighted <- diag(PCO.weighted$d[-K])
    U.unweighted <- PCO.unweighted$u[, -K] %*% sqrt(Lambdamat.unweighted)
    U.weighted <- PCO.weighted$u[, -K] %*% sqrt(Lambdamat.weighted)
    par(pty = "s")
    if (weight == "unweighted") {
        plot(U.unweighted[, 1:2], xlim = range(U.unweighted[, 
            1] * exp.factor), ylim = range(U.unweighted[, 2] * 
            exp.factor), asp = 1, pch = pch.means, col = pch.means.col, 
            cex = pch.means.size, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "")
        text(U.unweighted[, 1:2], labels = groupnames, pos = pos.labs, 
            cex = label.size)
        if (!is.null(X.new.samples)) {
            yvec.mat <- matrix(NA, nrow = nrow(X.new.samples), 
                ncol = K - 1)
            for (i in 1:nrow(X.new.samples)) {
                dvec <- dvec.fun(X.new.samples[i, ], X)
                dnplusone <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                  t(Gmat) %*% dvec
                yvec.mat[i, ] <- solve(Lambdamat.unweighted) %*% 
                  t(U.unweighted) %*% (dnplusone - apply(Deltabar, 
                  1, sum)/K)
            }
            points(yvec.mat[, 1:2], pch = pch.samples[samples.colvec], 
                col = pch.samples.col[samples.colvec], cex = pch.samples.size)
        }
    }
    if (weight == "weighted") {
        plot(U.weighted[, 1:2], xlim = range(U.weighted[, 1] * 
            exp.factor), ylim = range(U.weighted[, 2] * exp.factor), 
            asp = 1, pch = pch.means, col = pch.means.col, cex = pch.means.size, 
            xaxt = "n", yaxt = "n", xlab = "", ylab = "")
        text(U.weighted[, 1:2], labels = groupnames, pos = pos.labs, 
            cex = label.size)
    }
    list(T.ss = T.ss, B.ss = B.ss, W.ss = W.ss, nvec = nvec, 
        Dmat = Dmat, Dbar = Dbar, Deltabar = Deltabar)
}
