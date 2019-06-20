AOD.SS <-
function (X = Pine.data[, 2:8], class.vec = Pine.data[, 1], X.new.samples = scale(Pine.data[, 
    2:8]), scaled.mat = TRUE, e.vects = 1:2, dim.biplot = c(2, 
    1, 3), alpha = 0.95, ax = 1:ncol(X), ax.name.col = rep("black", 
    ncol(X)), ax.name.size = 0.75, ax.type = c("predictive", 
    "interpolative"), prediction.type = c("circle", "normal", 
    "back"), ax.col = list(ax.col = rep(8, ncol(X)), tickmarker.col = rep(8, 
    ncol(X)), marker.col = rep(1, ncol(X))), between = c(1, -1, 
    0, 1), between.columns = -1, cex.3d = 0.6, char.legend.size = c(1.2, 
    0.7), c.hull.n = 10, colour.scheme = NULL, colours = rep(1, 
    nrow(X)), columns = 1, dist = c("Pythagoras", "Clark", "SqrtL1"), 
    exp.factor = 1.2, expand.markervalsL = c(0, 0, 0.5, 0.2, 
        0.2, 0.2, 0.2), expand.markervalsR = c(1, 0, 0, 0, 0, 
        0.2, -0.155), label = TRUE, label.size = 0.8, large.scale = FALSE, 
    legend.type = c(means = FALSE, samples = TRUE, bags = FALSE), 
    line.length = c(1, 1), line.size = 2.5, line.type = 1, line.width = 1, 
    markers = TRUE, marker.size = 0.5, max.num = 2500, means.plot = TRUE, 
    n.int = c(10, 25, 10, 10, 10, 5, 5), num.points = 100, offset = rep(0.5, 
        4), ort.lty = 1, parlegendmar = c(3, 1, 3, 1), parplotmar = rep(3, 
        4), pch.means = rep(15, length(levels(class.vec))), pch.means.size = 1.5, 
    pch.means.col = UBcolours[1:5], pch.new = 1, pch.new.cols = "black", 
    pch.new.size = 1, pch.new.labels = NULL, pch.new.labels.size = 0.8, 
    pch.samples = rep(16, nrow(X)), pch.samples.size = 1, pch.samples.col = UBcolours[1:5], 
    plot.class.means = TRUE, pos = c("Orthog", "Hor", "Paral"), 
    pos.labs = 1, predictions.mean = NULL, predictions.sample = NULL, 
    quality.print = FALSE, specify.bags = 1:5, scale.axes = FALSE, 
    straight = "FALSE", specify.classes = 1:length(levels(class.vec)), 
    text.width.mult = 1, tra.width = 1, Title = "", Tukey.median = FALSE, 
    weight = c("unweighted", "weighted"), Wmat = "Default", zoomval = NULL, 
    epsilon = 0.001, select.origin = FALSE, ...) 
{
    G <- indmat(class.vec)
    dist <- dist[1]
    weight <- weight[1]
    prediction.type <- prediction.type[1]
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
    ynew.mat <- NULL
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
    D.fun <- function(X) {
        n <- nrow(X)
        D <- matrix(0, nrow = n, ncol = n)
        for (i in 1:(n - 1)) for (j in (i + 1):n) D[i, j] <- -0.5 * 
            sum(dist.expr(X[i, ], X[j, ]))
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
    dnplus1 <- function(X, x.new) {
        apply(X, 1, function(x, x.new) -0.5 * sum(dist.expr(x, 
            x.new)), x.new = x.new)
    }
    ddmu.expr <- function(x) {
        if (dist == "Pythagoras" || dist == "Clark") {
            dist.text <- deparse(dist.expr)[3]
            dist.text <- paste("-0.5", dist.text, sep = "*")
            dist.text <- gsub("xik", x, dist.text)
            dist.text <- gsub("xjk", "mu", dist.text)
            dist.expression <- parse(text = dist.text)
            AA <- D(dist.expression, "mu")
        }
        if (dist == "SqrtL1") 
            AA <- expression(-(sign(x - mu)))
        return(AA)
    }
    d2dmu2.expr <- function(x) {
        first <- ddmu.expr(x)
        D(first, "mu")
    }
    create.list.groupmat <- function(X, G) {
        indmat.logical <- apply(G, 2, function(x) as.logical(x))
        K <- ncol(G)
        list.groupmat <- vector("list", K)
        list.groupmat <- lapply(1:K, function(y) list.groupmat[[y]] <- X[indmat.logical[, 
            y], ])
        if (K == nrow(X)) 
            for (i in 1:K) list.groupmat[[i]] <- matrix(list.groupmat[[i]], 
                ncol = ncol(X))
        names(list.groupmat) <- colnames(G)
        list.groupmat
    }
    ax.type <- ax.type[1]
    pos <- pos[1]
    old.par <- par(no.readonly = TRUE)
    par(pty = "s", mar = parplotmar)
    X <- as.matrix(X)
    unscaled.X <- X
    means <- rep(0, ncol(X))
    sd <- rep(1, ncol(X))
    if (scaled.mat) 
        X <- scale(X)
    else {
        means <- rep(0, ncol(X))
        sd <- rep(1, ncol(X))
    }
    n <- nrow(X)
    p <- ncol(X)
    K <- ncol(G)
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste(1:n), paste("V", 1:p, sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste(1:n)
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
    groupnames <- colnames(G)
    samples.colvec <- apply(G, 1, function(x) (1:ncol(G))[x == 
        1])
    means.mat <- apply(X, 2, function(x) tapply(x, samples.colvec, 
        mean))
    Nmat <- t(G) %*% G
    nvec <- diag(Nmat)
    Dmat <- D.fun(X)
    Dbar <- solve(Nmat) %*% t(G) %*% Dmat %*% G %*% solve(Nmat)
    Deltamat <- Dbar - 0.5 * (diag(diag(Dbar)) %*% matrix(1, 
        nrow = K, ncol = K) + matrix(1, nrow = K, ncol = K) %*% 
        diag(diag(Dbar)))
    W.ss.alt <- -sum((apply(G, 2, function(x, D) {
        t(x) %*% D %*% x
    }, D = Dmat))/nvec)
    T.ss <- -sum(Dmat)/sum(nvec)
    B.ss <- -(t(nvec) %*% Deltamat %*% nvec)/sum(nvec)
    W.ss <- T.ss - B.ss
    centmat.n <- (diag(n) - matrix(1, ncol = n, nrow = n)/n)
    centmat.K <- (diag(K) - matrix(1, ncol = K, nrow = K)/K)
    centmat.nvec <- diag(K) - matrix(nvec, ncol = K, nrow = K, 
        byrow = TRUE)/sum(nvec)
    B.mat <- centmat.n %*% Dmat %*% centmat.n
    DoubCentDeltamat.unweighted <- centmat.K %*% Deltamat %*% 
        centmat.K
    DoubCentDeltamat.weighted <- centmat.nvec %*% Deltamat %*% 
        t(centmat.nvec)
    PCO.unweighted <- svd(DoubCentDeltamat.unweighted)
    PCO.weighted <- svd(DoubCentDeltamat.weighted)
    Lambdamat.unweighted <- zapsmall(diag(PCO.unweighted$d))
    Lambdamat.unweighted.inverse <- ifelse(Lambdamat.unweighted > 
        0, 1/Lambdamat.unweighted, 0)
    Lambdamat.weighted <- zapsmall(diag(PCO.weighted$d))
    Lambdamat.weighted.inverse <- ifelse(Lambdamat.weighted > 
        0, 1/Lambdamat.weighted, 0)
    U.unweighted <- PCO.unweighted$u %*% sqrt(Lambdamat.unweighted)
    if ((Wmat == "Default") | (Wmat == "default")) 
        Wmat <- Nmat
    U.weighted <- PCO.weighted$u %*% sqrt(Lambdamat.weighted)
    svd.1 <- svd(sqrt(Wmat) %*% U.weighted)
    Y.hat.2 <- U.weighted %*% svd.1$v
    Y.1.star <- centmat.nvec %*% U.unweighted
    Jmat <- diag(c(1, 1, rep(0, K - 2)))
    svd.3.Ybar2 <- svd(t(U.weighted) %*% Y.1.star)
    svd.3.Yhat2 <- svd(t(Y.hat.2) %*% Y.1.star)
    Qmat.Ybar2 <- svd.3.Ybar2$v %*% t(svd.3.Ybar2$u)
    Qmat.Yhat2 <- svd.3.Yhat2$v %*% t(svd.3.Yhat2$u)
    Y.1 <- U.unweighted
    Y.2 <- U.weighted
    Y.2.hat.Ybar2 <- Y.1.star %*% Qmat.Ybar2
    Y.2.hat.Yhat2 <- Y.1.star %*% Qmat.Yhat2
    plot(Y.hat.2[, e.vects], asp = 1)
    points(Y.2.hat.Ybar2[, e.vects], cex = 2)
    points(Y.2.hat.Yhat2[, e.vects], cex = 3, col = "red")
    Qmat <- Qmat.Yhat2
    U <- Y.1
    Y.samples.unweighted <- Lambdamat.unweighted.inverse %*% 
        t(U.unweighted) %*% sweep(t(G), 1, nvec, "/") %*% B.mat
    T.Delta.ss <- diag(round(Y.samples.unweighted %*% t(Y.samples.unweighted), 
        digits = 2))
    tempmat <- t(Y.samples.unweighted) %*% Qmat
    Y.samples.weighted <- t(tempmat)
    Pmat <- G %*% solve(t(G) %*% G) %*% t(G)
    T.Delta.ss.unweighted <- diag(round(Y.samples.unweighted %*% 
        t(Y.samples.unweighted), digits = 4))
    T.Delta.ss.weighted <- diag(round(Y.samples.weighted %*% 
        t(Y.samples.weighted), digits = 4))
    B.Delta.ss.unweighted <- diag(round(Y.samples.unweighted %*% 
        Pmat %*% t(Y.samples.unweighted), digits = 4))
    B.Delta.ss.weighted <- diag(round(Y.samples.weighted %*% 
        Pmat %*% t(Y.samples.weighted), digits = 4))
    W.Delta.ss.unweighted <- round(T.Delta.ss.unweighted - B.Delta.ss.unweighted, 
        digits = 4)
    W.Delta.ss.weighted <- round(T.Delta.ss.weighted - B.Delta.ss.weighted, 
        digits = 4)
    s.vec.unweighted <- (sum(Dmat)/(n^2)) * rep(1, n) - (2/n) * 
        apply(Dmat, 1, sum) - diag(t(Y.samples.unweighted) %*% 
        Y.samples.unweighted)
    s.vec.weighted <- (sum(Dmat)/(n^2)) * rep(1, n) - (2/n) * 
        apply(Dmat, 1, sum) - diag(t(Y.samples.weighted) %*% 
        Y.samples.weighted)
    W.OrthogDelta.ss.unweighted <- round(s.vec.unweighted, digits = 4)
    W.OrthogDelta.ss.weighted <- round(s.vec.weighted, digits = 4)
    W.OrthogDelta.totss.unweighted <- round(sum(s.vec.unweighted), 
        digits = 4)
    W.OrthogDelta.totss.weighted <- round(sum(s.vec.weighted), 
        digits = 4)
    W.OrthogDelta.indgroups.ss.unweighted <- round(t(G) %*% matrix(s.vec.unweighted, 
        ncol = 1), digits = 4)
    W.OrthogDelta.indgroups.ss.weighted <- round(t(G) %*% matrix(s.vec.weighted, 
        ncol = 1), digits = 4)
    out.ss <- list(T.ss = T.ss, B.ss = B.ss, W.ss = W.ss, T.Delta.ss.unweighted = T.Delta.ss.unweighted, 
        B.Delta.ss.unweighted = B.Delta.ss.unweighted, W.Delta.ss.unweighted = W.Delta.ss.unweighted, 
        W.OrthogDelta.ss.unweighted = W.OrthogDelta.ss.unweighted, 
        W.OrthogDelta.totss.unweighted = W.OrthogDelta.totss.unweighted, 
        W.OrthogDelta.indgroups.ss.unweighted = W.OrthogDelta.indgroups.ss.unweighted, 
        T.Delta.ss.weighted = T.Delta.ss.weighted, B.Delta.ss.weighted = B.Delta.ss.weighted, 
        W.Delta.ss.weighted = W.Delta.ss.weighted, W.OrthogDelta.ss.weighted = W.OrthogDelta.ss.weighted, 
        W.OrthogDelta.totss.weighted = W.OrthogDelta.totss.weighted, 
        W.OrthogDelta.indgroups.ss.weighted = W.OrthogDelta.indgroups.ss.weighted)
    if (weight == "unweighted") 
        eigval <- diag(t(Y.1) %*% Y.1)
    if (weight == "weighted") 
        eigval <- diag(t(Y.2) %*% Y.2)
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
        "%")
    if (weight == "unweighted") 
        Z <- Y.1[, e.vects]
    if (weight == "weighted") {
        Z <- Y.2.hat.Yhat2
        Z <- Z[, e.vects]
    }
    if (weight == "unweighted") {
        plot(U.unweighted[, e.vects], xlim = range(U.unweighted[, 
            e.vects[1]] * exp.factor), ylim = range(U.unweighted[, 
            e.vects[2]] * exp.factor), asp = 1, pch = pch.means, 
            col = pch.means.col, cex = pch.means.size, xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "")
        text(U.unweighted[, e.vects], labels = groupnames, pos = pos.labs, 
            cex = label.size)
        if (!is.null(X.new.samples)) {
            yvec.mat <- matrix(NA, nrow = nrow(X.new.samples), 
                ncol = K)
            for (i in 1:nrow(X.new.samples)) {
                dvec <- dvec.fun(X.new.samples[i, ], X)
                dnplusone <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                  t(G) %*% dvec
                yvec.mat[i, ] <- Lambdamat.unweighted.inverse %*% 
                  t(U.unweighted) %*% (dnplusone - apply(Deltamat, 
                  1, sum)/K)
            }
            points(yvec.mat[, e.vects, drop = FALSE], pch = pch.samples[samples.colvec], 
                col = pch.samples.col[samples.colvec], cex = pch.samples.size)
        }
    }
    if (weight == "weighted") {
        plot(Y.2.hat.Yhat2[, e.vects], xlim = range(Y.2.hat.Yhat2[, 
            e.vects[1]] * exp.factor), ylim = range(Y.2.hat.Yhat2[, 
            e.vects[2]] * exp.factor), asp = 1, pch = pch.means, 
            col = pch.means.col, cex = pch.means.size, xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "")
        text(Y.2.hat.Yhat2[, e.vects], labels = groupnames, pos = pos.labs, 
            cex = label.size)
        if (!is.null(X.new.samples)) {
            yvec.mat <- matrix(NA, nrow = nrow(X.new.samples), 
                ncol = K)
            for (i in 1:nrow(X.new.samples)) {
                dvec <- dvec.fun(X.new.samples[i, ], X)
                dnplusone <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                  t(G) %*% dvec
                yvec.mat[i, ] <- Lambdamat.unweighted.inverse %*% 
                  t(U.unweighted) %*% (dnplusone - apply(Deltamat, 
                  1, sum)/K)
            }
            yvec.mat <- (yvec.mat - matrix(nvec, nrow = nrow(X.new.samples), 
                ncol = K, byrow = TRUE) %*% Y.1/sum(nvec)) %*% 
                Qmat
            points(yvec.mat[, e.vects, drop = FALSE], pch = pch.samples[samples.colvec], 
                col = pch.samples.col[samples.colvec], cex = pch.samples.size)
        }
        points(x = 0, y = 0, pch = 3, cex = 3)
        origin.weight <- (rep(0, K) - matrix(nvec, nrow = 1, 
            ncol = K, byrow = TRUE) %*% Y.1/sum(nvec)) %*% Qmat
        points(x = origin.weight[1], y = origin.weight[2], pch = 4, 
            cex = 5, col = "blue")
    }
    means.s <- apply(unscaled.X, 2, mean)
    d.m <- D.fun(rbind(unscaled.X, 0))[n + 1, 1:n]
    s.vec <- rep(1/n, n)
    {
        dvec.m <- dvec.fun(means.s, X)
        dnplusone.m <- -0.5 * diag(Dbar) + solve(Nmat) %*% t(G) %*% 
            dvec.m
        O.co <- as.vector(Lambdamat.unweighted.inverse[e.vects, 
            e.vects] %*% t(U.unweighted[, e.vects]) %*% (dnplusone.m - 
            apply(Deltamat, 1, sum)/K))
    }
    lambda <- t(Y.1) %*% Y.1
    lambda.inv <- zapsmall(diag(ifelse(zapsmall(diag(lambda)) > 
        0, 1/diag(lambda), 0)))
    Z <- data.frame(Z, pch.sampl = pch.means, colr = as.character(pch.means.col), 
        line.type = line.type[1], stringsAsFactors = FALSE)
    dimnames(Z)[[1]] <- colnames(G)
    list.groupmat <- create.list.groupmat(unscaled.X, G)
    means.s <- apply(unscaled.X, 2, mean)
    sd.s <- sqrt(apply(unscaled.X, 2, var))
    if (ax.type == "predictive") {
        z.axes <- lapply(1:p, function(k, expandL, expandR, unscaled.X, 
            Xgrouplist, means, sd, n.int, D, Y, lambda, s, e.vecs, 
            prediction.type) {
            number.points <- num.points
            interval <- expand.interval(range(unscaled.X[, k]), 
                left = expandL[k], right = expandR[k])
            std.markers <- pretty(interval, n = n.int[k])
            interval <- (std.markers - means[k])/sd[k]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            number.points <- length(axis.vals)
            n <- nrow(unscaled.X)
            p <- ncol(unscaled.X)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            mu <- axis.vals
            m <- length(mu)
            K <- length(Xgrouplist)
            ddmu.dnplus1.list <- vector("list", K)
            for (i in 1:K) ddmu.dnplus1.list[[i]] <- t(sapply(Xgrouplist[[i]][, 
                k], function(x) eval(ddmu.expr(x)))/nrow(Xgrouplist[[i]]))
            cat(paste("Busy with axis ", k, "\n"))
            ddmu.dnplus1.mat <- matrix(NA, nrow = K, ncol = ncol(ddmu.dnplus1.list[[i]]))
            for (i in 1:K) ddmu.dnplus1.mat[i, ] <- apply(ddmu.dnplus1.list[[i]], 
                2, sum)
            A.mat <- t(ddmu.dnplus1.mat) %*% U[, e.vects] %*% 
                lambda.inv[e.vects, e.vects]
            m.star.vec <- apply(ddmu.dnplus1.mat, 2, sum)/(-K)
            a1.sq.plus.a2.sq <- apply(A.mat, 1, function(a) sum(a^2))
            LL.mat <- cbind(A.mat, m.star.vec, A.mat/sqrt(a1.sq.plus.a2.sq), 
                m.star.vec/sqrt(a1.sq.plus.a2.sq))
            L.a <- function(mu.val, Xlist, k) {
                K <- length(Xlist)
                out.list <- vector("list", K)
                for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                  k], function(x, mu) eval(ddmu.expr(x)), mu = mu.val)
                ddmu.dnplus1.vec <- rep(NA, K)
                for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                t(ddmu.dnplus1.vec) %*% U[, e.vects] %*% lambda.inv[e.vects, 
                  e.vects]
            }
            L.m <- function(mu.val, Xlist, k) {
                K <- length(Xlist)
                out.list <- vector("list", K)
                for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                  k], function(x, mu) eval(ddmu.expr(x)), mu = mu.val)
                ddmu.dnplus1.vec <- rep(NA, K)
                for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                m.star <- sum(ddmu.dnplus1.vec)/(-K)
                m.star/sum(L.a(mu.val, Xlist, k)^2)
            }
            if (prediction.type == "circle") {
                axis.points[, 1:2] <- cbind(LL.mat[, 4] * LL.mat[, 
                  6], LL.mat[, 5] * LL.mat[, 6])
            }
            if (prediction.type == "back") {
                mat <- sapply(mu, function(mm, unscaled.X, k) {
                  my.vek <- rep(0, ncol(unscaled.X))
                  my.vek[k] <- mm
                  dnplus1(unscaled.X, my.vek)
                }, unscaled.X = unscaled.X, k = k)
                xsi.k <- solve(lambda[e.vects, e.vects]) %*% 
                  t(Y[, e.vects]) %*% (mat - apply(D, 1, sum)/n)
                axis.points[, 1:2] <- t(sapply(1:m, function(j, 
                  xsi.k, LL.mat) {
                  l.vec <- LL.mat[j, 4:5]
                  m.mu <- LL.mat[j, 6]
                  xsi <- xsi.k[, j]
                  (diag(2) - l.vec %*% t(l.vec)) %*% xsi + m.mu * 
                    l.vec
                }, xsi.k = xsi.k, LL.mat = LL.mat))
            }
            if (prediction.type == "normal") {
                epsilon <- epsilon
                pos <- order(abs(LL.mat[, 6]))[1]
                if (LL.mat[pos, 6] < epsilon) 
                  mu.0 <- axis.vals[pos]
                else {
                  interval.begin <- pos - 1
                  if (interval.begin < 1) 
                    interval.begin <- 1
                  interval.einde <- pos + 1
                  if (interval.einde > number.points) 
                    interval.einde <- number.points
                  uit <- uniroot(L.m, axis.vals[c(interval.begin, 
                    interval.einde)], Xlist = list.groupmat, 
                    k)
                  mu.0 <- uit$root
                }
                skep.integrand.1 <- function(mu, Xlist, k) {
                  sapply(mu, function(mu, Xmatlist, k) {
                    K <- length(Xmatlist)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                      k], function(x, mu) eval(ddmu.expr(x)), 
                      mu = mu)
                    ddmu.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    a.vec <- t(ddmu.dnplus1.vec) %*% U[, e.vects] %*% 
                      lambda.inv[e.vects, e.vects]
                    wortel <- sqrt(sum(a.vec^2))
                    l1.mu <- a.vec[1]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-K)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xmatlist[[i]][, 
                      k], function(x, mu) eval(d2dmu2.expr(x)), 
                      mu = mu)
                    d2.dmu2.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) d2.dmu2.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-K)
                    ddmu.a1 <- t(U[, e.vects[1]]) %*% d2.dmu2.dnplus1.vec/lambda[e.vects[1], 
                      e.vects[1]]
                    ddmu.a2 <- t(U[, e.vects[2]]) %*% d2.dmu2.dnplus1.vec/lambda[e.vects[2], 
                      e.vects[2]]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l2mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[2] - (ddmu.a2 * 
                      wortel - ddmu.wortel * a.vec[2]) * mster.mu)/(wortel * 
                      a.vec[2]^2)
                    l1.mu * ddmu.mmu.l2mu
                  }, Xmatlist = Xlist, k = k)
                }
                f1.int <- lapply(mu, function(mu, Xgrouplist,k) try(integrate(skep.integrand.1, mu.0, mu, 
                       Xgrouplist, k)$value, silent=TRUE), Xgrouplist = Xgrouplist,  k = k)
				f1.int <- lapply(f1.int, function(x) if(is.null(attr(x,"condition"))) x <- x
                          else x <- attr(x,"condition"))
                f1.int <- unlist(lapply(f1.int, function(x)ifelse(is.numeric(x),x,NA)))
	   
				axis.points[, 1] <- LL.mat[, 5] * f1.int
				
                skep.integrand.2 <- function(mu, Xlist, k) {
                  sapply(mu, function(mu, Xmatlist, k) {
                    K <- length(Xmatlist)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                      k], function(x, mu) eval(ddmu.expr(x)), 
                      mu = mu)
                    ddmu.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    a.vec <- t(ddmu.dnplus1.vec) %*% U[, e.vects] %*% 
                      lambda.inv[e.vects, e.vects]
                    wortel <- sqrt(sum(a.vec^2))
                    l2.mu <- a.vec[2]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-K)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xmatlist[[i]][, 
                      k], function(x, mu) eval(d2dmu2.expr(x)), 
                      mu = mu)
                    d2.dmu2.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) d2.dmu2.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-K)
                    ddmu.a1 <- t(U[, e.vects[1]]) %*% d2.dmu2.dnplus1.vec/lambda[e.vects[1], 
                      e.vects[1]]
                    ddmu.a2 <- t(U[, e.vects[2]]) %*% d2.dmu2.dnplus1.vec/lambda[e.vects[2], 
                      e.vects[2]]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l1mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[1] - (ddmu.a1 * 
                      wortel - ddmu.wortel * a.vec[1]) * mster.mu)/(wortel * 
                      a.vec[1]^2)
                    l2.mu * ddmu.mmu.l1mu
                  }, Xmatlist = Xlist, k = k)
                }
				
				f2.int <- lapply(mu, function(mu, Xgrouplist,k) try(integrate(skep.integrand.2, mu.0, mu, 
                       Xgrouplist, k)$value, silent=TRUE), Xgrouplist = Xgrouplist,  k = k)
				f2.int <- lapply(f2.int, function(x) if(is.null(attr(x,"condition"))) x <- x
                          else x <- attr(x,"condition"))
                f2.int <- unlist(lapply(f2.int, function(x)ifelse(is.numeric(x),x,NA)))
				
                axis.points[, 2] <- LL.mat[, 4] * f2.int
            }
            axis.points[, 3] <- axis.vals * sd[k] + means[k]
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, expandL = expand.markervalsL, 
            expandR = expand.markervalsR, Xgrouplist = list.groupmat, 
            means = means, sd = sd, n.int = n.int, D = D, Y = U, 
            lambda = lambda, s = s, e.vecs = e.vects, prediction.type = prediction.type)
        z.axes <- sapply(z.axes, function(x) x <- na.omit(x))
    }
    yvec.mat <- NULL
    if (!is.null(X.new.samples)) {
        if (weight == "unweighted") {
            yvec.mat <- matrix(NA, nrow = nrow(X.new.samples), 
                ncol = K)
            for (i in 1:nrow(X.new.samples)) {
                dvec <- dvec.fun(X.new.samples[i, ], X)
                dnplusone <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                  t(G) %*% dvec
                yvec.mat[i, ] <- Lambdamat.unweighted.inverse %*% 
                  t(U.unweighted) %*% (dnplusone - apply(Deltamat, 
                  1, sum)/K)
            }
        }
        if (weight == "weighted") {
            plot(Y.2.hat.Yhat2[, e.vects], xlim = range(Y.2.hat.Yhat2[, 
                e.vects[1]] * exp.factor), ylim = range(Y.2.hat.Yhat2[, 
                e.vects[2]] * exp.factor), asp = 1, pch = pch.means, 
                col = pch.means.col, cex = pch.means.size, xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "")
            text(Y.2.hat.Yhat2[, e.vects], labels = groupnames, 
                pos = pos.labs, cex = label.size)
            if (!is.null(X.new.samples)) {
                yvec.mat <- matrix(NA, nrow = nrow(X.new.samples), 
                  ncol = K)
                for (i in 1:nrow(X.new.samples)) {
                  dvec <- dvec.fun(X.new.samples[i, ], X)
                  dnplusone <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                    t(G) %*% dvec
                  yvec.mat[i, ] <- Lambdamat.unweighted.inverse %*% 
                    t(U.unweighted) %*% (dnplusone - apply(Deltamat, 
                    1, sum)/K)
                }
                yvec.mat <- (yvec.mat - matrix(nvec, nrow = nrow(X.new.samples), 
                  ncol = K, byrow = TRUE) %*% Y.1/sum(nvec)) %*% 
                  Qmat
                points(yvec.mat[, e.vects, drop = FALSE], pch = pch.samples[samples.colvec], 
                  col = pch.samples.col[samples.colvec], cex = pch.samples.size)
            }
            points(x = 0, y = 0, pch = 3, cex = 3)
            origin.weight <- (rep(0, K) - matrix(nvec, nrow = 1, 
                ncol = K, byrow = TRUE) %*% Y.1/sum(nvec)) %*% 
                Qmat
            points(x = origin.weight[1], y = origin.weight[2], 
                pch = 4, cex = 4, col = "blue")
        }
    }
    if (ax.type == "interpolative") {
        z.axes <- lapply(1:p, function(j, unscaled.X, expandL, 
            expandR, means, means.2, sd, n.int, D, Y, lambda.inv, 
            dist.func, e.vecs, ax, K) {
            number.points <- num.points
            interval <- expand.interval(range(unscaled.X[, j]), 
                left = expandL[j], right = expandR[j])
            std.markers <- pretty(interval, n = n.int[j])
            interval <- (std.markers - means[j])/sd[j]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            number.points <- length(axis.vals)
            n <- nrow(unscaled.X)
            mat <- matrix(rep(means.2, rep(length(axis.vals), 
                ncol(unscaled.X))), ncol = ncol(unscaled.X), 
                nrow = length(axis.vals))
            mat[, j] <- axis.vals
            mat <- rbind(unscaled.X, mat)
            D.nuut <- dist.func(mat)
            dnplus1.mat <- D.nuut[1:n, -(1:n)]
            deltaK.mat <- -0.5 * (matrix(diag(Dbar), ncol = 1) %*% 
                matrix(1, nrow = 1, ncol = ncol(dnplus1.mat)) - 
                2 * solve(Nmat) %*% t(G) %*% dnplus1.mat)
            as.embed <- apply(deltaK.mat, 2, function(deltaK, 
                D, Y, lambda.inv, K) lambda.inv %*% t(Y) %*% 
                (deltaK - apply(D, 1, sum, D = D)/K), D = D, 
                Y = Y, K = K, lambda.inv = lambda.inv)
            axis.points <- cbind(t(as.embed)[, e.vecs], 0, 0)
            axis.points[, 3] <- axis.vals * sd[j] + means[j]
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            if (scale.axes) {
                dvec.m <- dvec.fun(means.s, X)
                dnplusone.m <- -0.5 * diag(Dbar) + solve(Nmat) %*% 
                  t(G) %*% dvec.m
                O.co <- as.vector(Lambdamat.unweighted.inverse[e.vecs, 
                  e.vecs] %*% t(U.unweighted[, e.vecs]) %*% (dnplusone.m - 
                  apply(Deltamat, 1, sum)/K))
                axis.points[, 1:2] <- scale(axis.points[, 1:2], 
                  center = O.co, scale = FALSE)
                axis.points[, 1:2] <- length(ax) * axis.points[, 
                  1:2]
                axis.points[, 1] <- axis.points[, 1] + O.co[1]
                axis.points[, 2] <- axis.points[, 2] + O.co[2]
            }
            return(axis.points)
        }, unscaled.X = unscaled.X, expandL = expand.markervalsL, 
            expandR = expand.markervalsR, means = means, means.2 = means.s, 
            sd = sd, n.int = n.int, D = Deltamat, Y = U, lambda.inv = lambda.inv, 
            dist.func = D.fun, e.vecs = e.vects, ax = ax, K = K)
    }
    temp.out <<- list(Z = Z, z.axes = z.axes, z.axes.names = colnames(X), 
        yvecnew.mat = yvec.mat, predictions.sample = predictions.sample)
    if (weight == "weighted") {
        for (i in 1:length(z.axes)) {
            zmat <- matrix(0, nrow = nrow(z.axes[[i]]), ncol = K)
            zmat[, 1:2] <- z.axes[[i]][, 1:2]
            z.axes[[i]][, 1:2] <- ((zmat - matrix(nvec, nrow = nrow(zmat), 
                ncol = K, byrow = TRUE) %*% Y.1/sum(nvec)) %*% 
                Qmat)[, e.vects]
        }
    }
    if (!is.null(X.new.samples)) {
        Z.means.mat <- Z
        Z <- data.frame(yvec.mat[, e.vects, drop = FALSE], pch = pch.samples[samples.colvec], 
            col = pch.samples.col[samples.colvec], cex = pch.samples.size)
    }
    else Z.means.mat <- NULL
    dev.new()
    par(pty = "s", mar = parplotmar)
    uit <- drawbipl.nonlin.AOD(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
        p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
        ax.name.col = ax.name.col, alpha = alpha, pch.means = pch.means, 
        pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
        specify.bags = specify.bags, label = label, markers = markers, 
        Title = Title, large.scale = large.scale, means.plot = means.plot, 
        specify.classes = specify.classes, Tukey.median = Tukey.median, 
        Z.means.mat = Z.means.mat, offset = offset, pos = pos, 
        strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
        marker.size = marker.size, label.size = label.size, exp.factor = exp.factor, 
        line.width = line.width, class.vec = class.vec, parplotmar = parplotmar, 
        predictions.sample = predictions.sample, predictions.mean = predictions.mean, 
        ort.lty = ort.lty, tra.width = tra.width, straight = straight, 
        ...)
    points(0, 0, pch = 3, cex = 5, col = "light grey")
    if (select.origin) {
        if (ax.type != "predictive") 
            stop("Select.origin only available for predictive biplot axes \n")
        cat("Move cursor where you want origin and press left mouse button \n")
        flush.console()
        origin.pos <- locator(1)
        dev.off()
        z.axes <- lapply(1:p, function(k, expandL, expandR, unscaled.X, 
            Xgrouplist, means, sd, n.int, D, Y, lambda, s, e.vecs, 
            prediction.type) {
            number.points <- num.points
            interval <- expand.interval(range(unscaled.X[, k]), 
                left = expandL[k], right = expandR[k])
            std.markers <- pretty(interval, n = n.int[k])
            interval <- (std.markers - means[k])/sd[k]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            number.points <- length(axis.vals)
            n <- nrow(unscaled.X)
            p <- ncol(unscaled.X)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            mu <- axis.vals
            m <- length(mu)
            K <- length(Xgrouplist)
            ddmu.dnplus1.list <- vector("list", K)
            for (i in 1:K) ddmu.dnplus1.list[[i]] <- t(sapply(Xgrouplist[[i]][, 
                k], function(x) eval(ddmu.expr(x)))/nrow(Xgrouplist[[i]]))
            cat(paste("Busy with axis ", k, "\n"))
            ddmu.dnplus1.mat <- matrix(NA, nrow = K, ncol = ncol(ddmu.dnplus1.list[[i]]))
            for (i in 1:K) ddmu.dnplus1.mat[i, ] <- apply(ddmu.dnplus1.list[[i]], 
                2, sum)
            A.mat <- t(ddmu.dnplus1.mat) %*% U[, e.vecs] %*% 
                lambda.inv[e.vecs, e.vecs]
            m.star.vec <- apply(ddmu.dnplus1.mat, 2, sum)/(-K)
            a1.sq.plus.a2.sq <- apply(A.mat, 1, function(a) sum(a^2))
            LL.mat <- cbind(A.mat, m.star.vec, A.mat/sqrt(a1.sq.plus.a2.sq), 
                m.star.vec/sqrt(a1.sq.plus.a2.sq))
            L.a <- function(mu.val, Xlist, k) {
                K <- length(Xlist)
                out.list <- vector("list", K)
                for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                  k], function(x, mu) eval(ddmu.expr(x)), mu = mu.val)
                ddmu.dnplus1.vec <- rep(NA, K)
                for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                t(ddmu.dnplus1.vec) %*% U[, e.vecs] %*% lambda.inv[e.vecs, 
                  e.vecs]
            }
            L.m <- function(mu.val, Xlist, k) {
                K <- length(Xlist)
                out.list <- vector("list", K)
                for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                  k], function(x, mu) eval(ddmu.expr(x)), mu = mu.val)
                ddmu.dnplus1.vec <- rep(NA, K)
                for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                m.star <- sum(ddmu.dnplus1.vec)/(-K)
                m.star/sum(L.a(mu.val, Xlist, k)^2)
            }
            if (prediction.type == "circle") {
                axis.points[, 1:2] <- cbind(LL.mat[, 4] * LL.mat[, 
                  6], LL.mat[, 5] * LL.mat[, 6])
            }
            if (prediction.type == "back") {
                mat <- sapply(mu, function(mm, unscaled.X, k) {
                  my.vek <- rep(0, ncol(unscaled.X))
                  my.vek[k] <- mm
                  dnplus1(unscaled.X, my.vek)
                }, unscaled.X = unscaled.X, k = k)
                xsi.k <- solve(lambda[e.vecs, e.vecs]) %*% t(Y[, 
                  e.vecs]) %*% (mat - apply(D, 1, sum)/n)
                axis.points[, 1:2] <- t(sapply(1:m, function(j, 
                  xsi.k, LL.mat) {
                  l.vec <- LL.mat[j, 4:5]
                  m.mu <- LL.mat[j, 6]
                  xsi <- xsi.k[, j]
                  (diag(2) - l.vec %*% t(l.vec)) %*% xsi + m.mu * 
                    l.vec
                }, xsi.k = xsi.k, LL.mat = LL.mat))
            }
            if (prediction.type == "normal") {
                epsilon <- epsilon
                pos <- order(abs(LL.mat[, 6]))[1]
                if (LL.mat[pos, 6] < epsilon) 
                  mu.0 <- axis.vals[pos]
                else {
                  interval.begin <- pos - 1
                  if (interval.begin < 1) 
                    interval.begin <- 1
                  interval.einde <- pos + 1
                  if (interval.einde > number.points) 
                    interval.einde <- number.points
                  uit <- uniroot(L.m, axis.vals[c(interval.begin, 
                    interval.einde)], Xlist = list.groupmat, 
                    k)
                  mu.0 <- uit$root
                }
                skep.integrand.1 <- function(mu, Xlist, k) {
                  sapply(mu, function(mu, Xmatlist, k) {
                    K <- length(Xmatlist)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                      k], function(x, mu) eval(ddmu.expr(x)), 
                      mu = mu)
                    ddmu.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    a.vec <- t(ddmu.dnplus1.vec) %*% U[, e.vecs] %*% 
                      lambda.inv[e.vecs, e.vecs]
                    wortel <- sqrt(sum(a.vec^2))
                    l1.mu <- a.vec[1]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-K)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xmatlist[[i]][, 
                      k], function(x, mu) eval(d2dmu2.expr(x)), 
                      mu = mu)
                    d2.dmu2.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) d2.dmu2.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-K)
                    ddmu.a1 <- t(U[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(U[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l2mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[2] - (ddmu.a2 * 
                      wortel - ddmu.wortel * a.vec[2]) * mster.mu)/(wortel * 
                      a.vec[2]^2)
                    l1.mu * ddmu.mmu.l2mu
                  }, Xmatlist = Xlist, k = k)
                }
          				  
				f1.int <- lapply(mu, function(mu, Xgrouplist,k) try(integrate(skep.integrand.1, mu.0, mu, 
                       Xgrouplist, k)$value, silent=TRUE), Xgrouplist = Xgrouplist,  k = k)
				f1.int <- lapply(f1.int, function(x) if(is.null(attr(x,"condition"))) x <- x
                          else x <- attr(x,"condition"))
                f1.int <- unlist(lapply(f1.int, function(x)ifelse(is.numeric(x),x,NA)))  
				  
                axis.points[, 1] <- LL.mat[, 5] * f1.int
                skep.integrand.2 <- function(mu, Xlist, k) {
                  sapply(mu, function(mu, Xmatlist, k) {
                    K <- length(Xmatlist)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xlist[[i]][, 
                      k], function(x, mu) eval(ddmu.expr(x)), 
                      mu = mu)
                    ddmu.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) ddmu.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    a.vec <- t(ddmu.dnplus1.vec) %*% U[, e.vecs] %*% 
                      lambda.inv[e.vecs, e.vecs]
                    wortel <- sqrt(sum(a.vec^2))
                    l2.mu <- a.vec[2]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-K)
                    out.list <- vector("list", K)
                    for (i in 1:K) out.list[[i]] <- sapply(Xmatlist[[i]][, 
                      k], function(x, mu) eval(d2dmu2.expr(x)), 
                      mu = mu)
                    d2.dmu2.dnplus1.vec <- rep(NA, K)
                    for (i in 1:K) d2.dmu2.dnplus1.vec[i] <- sum(out.list[[i]])/nrow(Xlist[[i]])
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-K)
                    ddmu.a1 <- t(U[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(U[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l1mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[1] - (ddmu.a1 * 
                      wortel - ddmu.wortel * a.vec[1]) * mster.mu)/(wortel * 
                      a.vec[1]^2)
                    l2.mu * ddmu.mmu.l1mu
                  }, Xmatlist = Xlist, k = k)
                }
				  
				f2.int <- lapply(mu, function(mu, Xgrouplist,k) try(integrate(skep.integrand.2, mu.0, mu, 
                       Xgrouplist, k)$value, silent=TRUE), Xgrouplist = Xgrouplist,  k = k)
				f2.int <- lapply(f2.int, function(x) if(is.null(attr(x,"condition"))) x <- x
                          else x <- attr(x,"condition"))
                f2.int <- unlist(lapply(f2.int, function(x)ifelse(is.numeric(x),x,NA))) 
				
                axis.points[, 2] <- LL.mat[, 4] * f2.int
            }
            axis.points[, 3] <- axis.vals * sd[k] + means[k]
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, expandL = expand.markervalsL, 
            expandR = expand.markervalsR, Xgrouplist = list.groupmat, 
            means = means, sd = sd, n.int = n.int, D = D, Y = U, 
            lambda = lambda, s = s, e.vecs = e.vects, prediction.type = prediction.type)
        z.axes <- sapply(z.axes, function(x) x <- na.omit(x))
        par(pty = "s", mar = parplotmar)
        uit <- drawbipl.nonlin.AOD(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
            ax.name.col = ax.name.col, alpha = alpha, pch.means = pch.means, 
            pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
            specify.bags = specify.bags, label = label, markers = markers, 
            Title = Title, large.scale = large.scale, means.plot = means.plot, 
            specify.classes = specify.classes, Tukey.median = Tukey.median, 
            Z.means.mat = Z.means.mat, offset = offset, pos = pos, 
            strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, parplotmar = parplotmar, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, ort.lty = ort.lty, 
            straight = straight, tra.width = tra.width, ...)
        points(0, 0, pch = 3, cex = 5, col = "light grey")
    }
    if (!is.null(zoomval)) {
        zoomval <- zoom(zoomval)
        drawbipl.nonlin.AOD(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
            ax.name.col = ax.name.col, alpha = alpha, pch.means = pch.means, 
            pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
            specify.bags = specify.bags, label = label, markers = markers, 
            Title = Title, large.scale = large.scale, means.plot = means.plot, 
            specify.classes = specify.classes, Tukey.median = Tukey.median, 
            Z.means.mat = Z.means.mat, offset = offset, pos = pos, 
            strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, parplotmar = parplotmar, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, ort.lty = ort.lty, 
            straight = straight, tra.width = tra.width, zoomval = zoomval)
        points(0, 0, pch = 3, cex = 5, col = "light grey")
        if (!is.null(X.new.samples)) 
            points(yvec.mat[, 1:2, drop = FALSE], pch = pch.samples[samples.colvec], 
                col = pch.samples.col[samples.colvec], cex = pch.samples.size)
    }
    if (any(legend.type)) {
        dev.new()
        blegend.colchar(quality.print = quality.print, QualityOfDisplay = display.quality, 
            classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:K]), 
            pch.samples = as.vector(pch.samples[1:K]), colours = as.vector(pch.samples.col[1:K]), 
            colours.means = as.vector(pch.means.col[1:K]), line.type = as.vector(line.type[1:K]), 
            line.width = line.width, pch.samples.size = char.legend.size, 
            legend.type = legend.type, between = between, columns = columns, 
            parlegendmar = parlegendmar, line.size = line.size, 
            between.columns = between.columns, text.width.mult = text.width.mult)
    }
    list(D = Dmat, Z = Z, Z.axes = z.axes, e.vals = NULL, predictions = uit$predictions, 
        O.co = O.co, ynew.mat = yvec.mat, U = U, Qmat = Qmat, 
        quality = fit.quality, out.ss, axes.names = dimnames(X)[[2]])
}
