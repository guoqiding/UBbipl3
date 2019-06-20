CVA.predictivities <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    weightedCVA = c("weighted", "unweightedI", "unweightedCent"), 
    X.new.samples = NULL, X.new.vars = NULL) 
{
    weightedCVA <- weightedCVA[1]
    OverAllMeans <- apply(X, 2, mean)
    Axis.predictivities.new.1 <- NULL
    Axis.predictivities.new.2 <- NULL
    new.samples.pred <- NULL
    X <- as.matrix(X)
    unscaled.X <- X
    X.cent <- scale(X, center = TRUE, scale = FALSE)
    OverAllMeans <- apply(X, 2, mean)
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    K <- min(p, J - 1)
    Gmat <- G
    Nmat <- t(Gmat) %*% Gmat
    XBar.groups <- solve(Nmat) %*% t(Gmat) %*% X
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
    canon.group.means <- XcentBar.groups %*% Mmat
    lambdamat.2sided <- svd.2sided$Lambda.mat
    Uit.list.Axis.pred <- vector("list", p)
    Uit.list.Class.pred <- vector("list", p)
    Uit.vec.CVA.quality.origvar <- rep(NA, p)
    Uit.list.Axis.adeq <- vector("list", p)
    for (i in 1:p) {
        vec.temp <- rep(0, p)
        vec.temp[1:i] <- 1
        Jmat <- diag(vec.temp)
        XcentBarLVJ <- XcentBar.groups %*% Mmat %*% Jmat
        XcentLVJ <- X.cent %*% Mmat %*% Jmat
        XcentBarHat <- XcentBar.groups %*% Mmat %*% Jmat %*% 
            solve(Mmat)
        Uit.list.Axis.adeq[[i]] <- diag(diag(Mmat %*% Jmat %*% 
            t(Mmat)))/diag(diag(Mmat %*% t(Mmat)))
        Uit.list.Axis.pred[[i]] <- diag(diag(t(XcentBarHat) %*% 
            Cmat %*% XcentBarHat))/diag(diag(t(XcentBar.groups) %*% 
            Cmat %*% XcentBar.groups))
        Uit.list.Class.pred[[i]] <- diag(diag(Csqrt %*% XcentBarHat %*% 
            solve(Wmat) %*% t(XcentBarHat) %*% Csqrt))/diag(diag(Csqrt %*% 
            XcentBar.groups %*% solve(Wmat) %*% t(XcentBar.groups) %*% 
            Csqrt))
        Uit.vec.CVA.quality.origvar[i] <- 100 * sum(diag(t(XcentBarHat) %*% 
            Cmat %*% XcentBarHat))/sum(diag(t(XcentBar.groups) %*% 
            Cmat %*% XcentBar.groups))
    }
    names(Uit.vec.CVA.quality.origvar) <- paste("Dim", 1:p, sep = "_")
    Axis.Adequacies <- round(diag(Uit.list.Axis.adeq[[1]]), digits = 4)
    for (i in 2:p) Axis.Adequacies <- rbind(Axis.Adequacies, 
        round(diag(Uit.list.Axis.adeq[[i]]), digits = 4))
    dimnames(Axis.Adequacies) <- list(paste("Dim", 1:p, sep = "_"), 
        dimnames(X)[[2]])
    Axis.Predictivities <- round(diag(Uit.list.Axis.pred[[1]]), 
        digits = 4)
    for (i in 2:p) Axis.Predictivities <- rbind(Axis.Predictivities, 
        round(diag(Uit.list.Axis.pred[[i]]), digits = 4))
    dimnames(Axis.Predictivities) <- list(paste("Dim", 1:p, sep = "_"), 
        dimnames(X)[[2]])
    Class.Predictivities <- round(diag(Uit.list.Class.pred[[1]]), 
        digits = 4)
    for (i in 2:p) Class.Predictivities <- rbind(Class.Predictivities, 
        round(diag(Uit.list.Class.pred[[i]]), digits = 4))
    dimnames(Class.Predictivities) <- list(paste("Dim", 1:p, 
        sep = "_"), dimnames(G)[[2]])
    X.within <- (diag(n) - Gmat %*% (solve(Nmat)) %*% t(Gmat)) %*% 
        X.cent
    Uit.list.WithinGroup.pred <- vector("list", p)
    for (i in 1:p) {
        vec.temp <- rep(0, p)
        vec.temp[1:i] <- 1
        Jmat <- diag(vec.temp)
        X.within.hat <- X.within %*% Mmat %*% Jmat %*% solve(Mmat)
        Uit.list.WithinGroup.pred[[i]] <- diag(diag(t(X.within.hat) %*% 
            X.within.hat))/(diag(diag(t(X.within) %*% X.within)))
    }
    WithinGroup.axis.Predictivities <- round(diag(Uit.list.WithinGroup.pred[[1]]), 
        digits = 4)
    for (i in 2:p) WithinGroup.axis.Predictivities <- rbind(WithinGroup.axis.Predictivities, 
        round(diag(Uit.list.WithinGroup.pred[[i]]), digits = 4))
    dimnames(WithinGroup.axis.Predictivities) <- list(paste("Dim", 
        1:p, sep = "_"), dimnames(X)[[2]])
    Uit.list.WithinGroup.Sample.pred <- vector("list", p)
    for (i in 1:p) {
        vec.temp <- rep(0, p)
        vec.temp[1:i] <- 1
        Jmat <- diag(vec.temp)
        X.within.hat <- X.within %*% Mmat %*% Jmat %*% solve(Mmat)
        Uit.list.WithinGroup.Sample.pred[[i]] <- diag(diag(((X.within.hat) %*% 
            solve(Wmat) %*% t(X.within.hat))))/diag(((diag((X.within) %*% 
            solve(Wmat) %*% t(X.within)))))
    }
    WithinGroup.Sample.Predictivities <- round(diag(Uit.list.WithinGroup.Sample.pred[[1]]), 
        digits = 4)
    for (i in 2:p) WithinGroup.Sample.Predictivities <- rbind(WithinGroup.Sample.Predictivities, 
        round(diag(Uit.list.WithinGroup.Sample.pred[[i]]), digits = 4))
    dimnames(WithinGroup.Sample.Predictivities) <- list(paste("Dim", 
        1:p, sep = "_"), dimnames(X)[[1]])
    Quality.canvar <- (100 * cumsum(diag(lambdamat)))/sum(diag(lambdamat))
    names(Quality.canvar) <- paste("Dim", 1:length(Quality.canvar))
    if (!is.null(X.new.samples)) {
        if (!is.matrix(X.new.samples)) 
            stop("X.new.samples must be a t x p matrix \n")
        new.samples.pred <- vector("list", nrow(X.new.samples))
        x.new.star <- scale(X.new.samples, center = OverAllMeans, 
            scale = FALSE)
        for (aa in 1:nrow(X.new.samples)) {
            mat.temp <- matrix(NA, nrow = p, ncol = J)
            for (i in 1:p) {
                vec.temp <- rep(0, p)
                vec.temp[1:i] <- 1
                Jmat <- diag(vec.temp)
                new.samples.pred.row <- rep(NA, J)
                for (j in 1:J) {
                  group.mean.row <- XcentBar.groups[j, ]
                  x.tilde <- matrix(x.new.star[aa, ] - group.mean.row, 
                    nrow = p, ncol = 1)
                  new.samples.pred.row[j] <- (t(x.tilde) %*% 
                    Mmat %*% Jmat %*% t(Mmat) %*% x.tilde)/(t(x.tilde) %*% 
                    Mmat %*% t(Mmat) %*% x.tilde)
                }
                mat.temp[i, ] <- round(new.samples.pred.row, 
                  digits = 4)
            }
            dimnames(mat.temp) <- list(paste("Dim_", 1:p, sep = ""), 
                rownames(XcentBar.groups))
            new.samples.pred[[aa]] <- mat.temp
        }
    }
    if (!is.null(X.new.vars)) {
        if (!is.matrix(X.new.vars)) 
            stop("If X.new.vars is not NULL it must be in the form of a matrix \n")
        NewVars.means <- apply(X.new.vars, 2, mean)
        NewVars.cent <- scale(X.new.vars, center = TRUE, scale = FALSE)
        NewVarsBar.groups <- solve(Nmat) %*% t(Gmat) %*% NewVars.cent
        LambdaMinOne <- ifelse(lambdamat < 1e-10, 0, 1/lambdamat)
        b.regres <- LambdaMinOne %*% t(Mmat) %*% t(XcentBar.groups) %*% 
            Cmat %*% NewVarsBar.groups
        p <- ncol(X)
        Out.list.vars.new.1 <- vector("list", p)
        Out.list.vars.new.2 <- vector("list", p)
        for (i in 1:p) {
            Jr <- diag(c(rep(1, i), rep(0, p - i)))
            numer.mat <- t(b.regres) %*% Jr %*% lambdamat %*% 
                Jr %*% b.regres
            denom.mat.1 <- t(NewVarsBar.groups) %*% Cmat %*% 
                NewVarsBar.groups
            denom.mat.2 <- t(b.regres) %*% lambdamat %*% b.regres
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
    list(OverAllMeans = round(OverAllMeans, digits = 4), XBar.groups = round(XBar.groups, 
        digits = 4), XcentBar.groups = round(XcentBar.groups, 
        digits = 4), canon.group.means = canon.group.means, CVA.quality.canvar = round(Quality.canvar, 
        digits = 4), CVA.quality.origvar = round(Uit.vec.CVA.quality.origvar, 
        digits = 4), Axis.Adequacies = Axis.Adequacies, Axis.Predictivities = Axis.Predictivities, 
        Class.Predictivities = Class.Predictivities, WithinGroup.axis.Predictivities = WithinGroup.axis.Predictivities, 
        WithinGroup.Sample.Predictivities = WithinGroup.Sample.Predictivities, 
        new.samples.pred = new.samples.pred, Axis.predictivities.new.vars.1 = Axis.predictivities.new.1, 
        Axis.predictivities.new.vars.2 = Axis.predictivities.new.2)
}
