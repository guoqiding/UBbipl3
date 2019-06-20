PCA.predictivities <-
function (X, X.new.samples = NULL, X.new.vars = NULL, scaled.mat = FALSE) 
{
    if (!is.null(X.new.samples)) 
        if (!is.matrix(X.new.samples) | !identical(ncol(X), ncol(X.new.samples))) 
            stop("X.new.samples must be a matrix of size k x p \n")
    if (!is.null(X.new.vars)) 
        if (!is.matrix(X.new.vars) | !identical(nrow(X), nrow(X.new.vars))) 
            stop("X.new.vars must be a matrix of size n x k \n")
    means.original <- apply(X, 2, mean)
    if (identical(scaled.mat, TRUE)) {
        sds.original <- sqrt(apply(X, 2, var))
        X <- scale(X)
        if (!is.null(X.new.samples)) 
            X.new.samples <- scale(X.new.samples, scale = sds.original, 
                center = means.original)
        if (!is.null(X.new.vars)) 
            X.new.vars <- scale(X.new.vars, scale = TRUE, center = TRUE)
    }
    else {
        X <- scale(X, scale = FALSE, center = TRUE)
        if (!is.null(X.new.samples)) 
            X.new.samples <- scale(X.new.samples, scale = FALSE, 
                center = means.original)
        if (!is.null(X.new.vars)) 
            X.new.vars <- scale(X.new.vars, scale = FALSE, center = TRUE)
    }
    out.X <- svd(X)
    Sigma2 <- diag(out.X$d^2)
    Quality <- (100 * cumsum(diag(Sigma2)))/sum(diag(Sigma2))
    names(Quality) <- paste("Dim", 1:length(Quality), sep = "_")
    out.XX <- svd(X %*% t(X))
    Umat <- out.X$u
    Vmat <- out.X$v
    Umat.full <- out.XX$u
    Weights <- diag((Vmat %*% Sigma2 %*% t(Vmat)))/sum(Sigma2)
    p <- ncol(X)
    Out.list.samples <- vector("list", p)
    for (i in 1:p) {
        Jr <- diag(c(rep(1, i), rep(0, p - i)))
        numer.mat <- X %*% Vmat %*% Jr %*% t(Vmat) %*% t(X)
        denom.mat <- X %*% t(X)
        Out.list.samples[[i]] <- diag(diag(numer.mat)/diag(denom.mat))
    }
    Sample.predictivities.original <- round(diag(Out.list.samples[[1]]), 
        digits = 4)
    for (i in 2:p) Sample.predictivities.original <- rbind(Sample.predictivities.original, 
        round(diag(Out.list.samples[[i]]), digits = 4))
    dimnames(Sample.predictivities.original) <- list(paste("Dim", 
        1:p, sep = "_"), dimnames(X)[[1]])
    if (!is.null(X.new.samples)) {
        p <- ncol(X)
        Out.list.samples.new <- vector("list", p)
        for (i in 1:p) {
            Jr <- diag(c(rep(1, i), rep(0, p - i)))
            numer.mat <- X.new.samples %*% Vmat %*% Jr %*% t(Vmat) %*% 
                t(X.new.samples)
            denom.mat <- X.new.samples %*% t(X.new.samples)
            if (nrow(denom.mat) < 2) 
                Out.list.samples.new[[i]] <- numer.mat/denom.mat
            else Out.list.samples.new[[i]] <- diag(diag(numer.mat)/diag(denom.mat))
        }
        Sample.predictivities.new <- round(diag(Out.list.samples.new[[1]]), 
            digits = 4)
        for (i in 2:p) Sample.predictivities.new <- rbind(Sample.predictivities.new, 
            round(diag(Out.list.samples.new[[i]]), digits = 4))
        dimnames(Sample.predictivities.new) <- list(paste("Dim", 
            1:p, sep = "_"), dimnames(X.new.samples)[[1]])
    }
    p <- ncol(X)
    Out.list.vars <- vector("list", p)
    Out.list.Adeq <- vector("list", p)
    for (i in 1:p) {
        Jr <- diag(c(rep(1, i), rep(0, p - i)))
        numer.mat <- t(X) %*% Umat %*% Jr %*% t(Umat) %*% X
        denom.mat <- t(X) %*% X
        Out.list.vars[[i]] <- diag(diag(numer.mat)/diag(denom.mat))
        Out.list.Adeq[[i]] <- diag(diag(Vmat %*% Jr %*% t(Vmat)))
    }
    Axis.predictivities.original <- round(diag(Out.list.vars[[1]]), 
        digits = 4)
    for (i in 2:p) Axis.predictivities.original <- rbind(Axis.predictivities.original, 
        round(diag(Out.list.vars[[i]]), digits = 4))
    dimnames(Axis.predictivities.original) <- list(paste("Dim", 
        1:p, sep = "_"), dimnames(X)[[2]])
    Adequacies <- round(diag(Out.list.Adeq[[1]]), digits = 4)
    for (i in 2:p) Adequacies <- rbind(Adequacies, round(diag(Out.list.Adeq[[i]]), 
        digits = 4))
    dimnames(Adequacies) <- list(paste("Dim", 1:p, sep = "_"), 
        dimnames(X)[[2]])
    if (!is.null(X.new.vars)) {
        p <- ncol(X)
        Out.list.vars.new.1 <- vector("list", p)
        Out.list.vars.new.2 <- vector("list", p)
        for (i in 1:p) {
            Jr <- diag(c(rep(1, i), rep(0, p - i)))
            numer.mat <- t(X.new.vars) %*% Umat %*% Jr %*% t(Umat) %*% 
                X.new.vars
            denom.mat.1 <- t(X.new.vars) %*% X.new.vars
            denom.mat.2 <- t(X.new.vars) %*% Umat %*% t(Umat) %*% 
                X.new.vars
            if (nrow(denom.mat.1) < 2) {
                Out.list.vars.new.1[[i]] <- numer.mat/denom.mat.1
                Out.list.vars.new.2[[i]] <- numer.mat/denom.mat.2
            }
            else {
                Out.list.vars.new.1[[i]] <- diag(diag(numer.mat)/diag(denom.mat.1))
                Out.list.vars.new.2[[i]] <- diag(diag(numer.mat)/diag(denom.mat.2))
            }
        }
        Axis.predictivities.new.1 <- round(diag(Out.list.vars.new.1[[1]]), 
            digits = 4)
        for (i in 2:p) Axis.predictivities.new.1 <- rbind(Axis.predictivities.new.1, 
            round(diag(Out.list.vars.new.1[[i]]), digits = 4))
        dimnames(Axis.predictivities.new.1) <- list(paste("Dim", 
            1:p, sep = "_"), dimnames(X.new.vars)[[2]])
        Axis.predictivities.new.2 <- round(diag(Out.list.vars.new.2[[1]]), 
            digits = 4)
        for (i in 2:p) Axis.predictivities.new.2 <- rbind(Axis.predictivities.new.2, 
            round(diag(Out.list.vars.new.2[[i]]), digits = 4))
        dimnames(Axis.predictivities.new.2) <- list(paste("Dim", 
            1:p, sep = "_"), dimnames(X.new.vars)[[2]])
    }
    if (is.null(X.new.samples)) 
        Sample.predictivities.new <- NULL
    if (is.null(X.new.vars)) {
        Axis.predictivities.new.1 <- NULL
        Axis.predictivities.new.2 <- NULL
    }
    list(Quality = Quality, Weights = Weights, Sample.predictivities.original = Sample.predictivities.original, 
        Adequacies = Adequacies, Axis.predictivities.original = Axis.predictivities.original, 
        Sample.predictivities.new = Sample.predictivities.new, 
        Axis.predictivities.new.1 = Axis.predictivities.new.1, 
        Axis.predictivities.new.2 = Axis.predictivities.new.2)
}
