Nonlinbipl <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    X.new.samples = NULL, scaled.mat = FALSE, e.vects = 1:2, 
    dim.biplot = c(2, 1, 3), alpha = 0.95, ax = 1:ncol(X), ax.name.col = rep("black", 
        ncol(X)), ax.name.size = 0.75, ax.type = c("predictive", 
        "interpolative"), prediction.type = c("circle", "normal", 
        "back"), ax.col = list(ax.col = rep(8, ncol(X)), tickmarker.col = rep(8, 
        ncol(X)), marker.col = rep(1, ncol(X))), between = c(1, 
        -1, 0, 1), between.columns = -1, cex.3d = 0.6, char.legend.size = c(1.2, 
        0.7), c.hull.n = 10, colour.scheme = NULL, colours = rep(1, 
        nrow(X)), columns = 1, exp.factor = 1.2, label = TRUE, 
    label.size = 0.8, large.scale = FALSE, legend.type = c(means = FALSE, 
        samples = TRUE, bags = FALSE), line.length = c(1, 1), 
    line.size = 2.5, line.type = 1, line.width = 1, markers = TRUE, 
    marker.size = 0.5, max.num = 2500, means.plot = TRUE, n.int = rep(5, 
        ncol(X)), offset = rep(0.5, 4), ort.lty = 1, parlegendmar = c(3, 
        1, 3, 1), parplotmar = rep(3, 4), plot.class.means = TRUE, 
    pch.means = 0:10, pch.means.size = 1, pch.new = 1, pch.new.cols = "black", 
    pch.new.size = 1, pch.new.labels = NULL, pch.new.labels.size = 0.8, 
    pch.samples = rep(15, nrow(X)), pch.samples.size = 1, pos = c("Orthog", 
        "Hor", "Paral"), predictions.mean = NULL, predictions.sample = NULL, 
    quality.print = FALSE, specify.bags = dimnames(G)[[2]], specify.classes = dimnames(G)[[2]], 
    text.width.mult = 1, Title = "", Tukey.median = FALSE, dist = c("Pythagoras", 
        "Clark", "SqrtL1"), scale.axes = FALSE, straight = "FALSE", 
    num.points = 25, zoomval = NULL) 
{
    dist <- dist[1]
    prediction.type <- prediction.type[1]
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
    D.mat <- function(X) {
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
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste(1:n), paste("V", 1:p, sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste(1:n)
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
    s <- rep(1/n, n)
    one <- rep(1, n)
    N <- one %*% t(s)
    I <- diag(n)
    D <- D.mat(X)
    B <- (I - N) %*% D %*% (I - N)
    swd <- svd(B)
    U <- (swd$u %*% diag(swd$d^0.5))[, -n]
    lambda <- t(U) %*% U
    lambda.inv <- diag(ifelse(zapsmall(diag(lambda)) > 0, 1/diag(lambda), 
        0))
    eigval <- diag(lambda)
    eigval.r <- eigval[e.vects[1:2]]
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
        "%")
    means.s <- apply(unscaled.X, 2, mean)
    d.m <- D.mat(rbind(unscaled.X, 0))[n + 1, 1:n]
    O.co <- as.vector(lambda.inv[1:2, 1:2] %*% t(U[, 1:2]) %*% 
        (d.m - D %*% s))
    Z <- U[, e.vects[1:2]]
    Z <- data.frame(Z, pch.sampl = pch.samples[1], colr = as.character(colours), 
        line.type = line.type[1], stringsAsFactors = FALSE)
    dimnames(Z)[[1]] <- dimnames(X)[[1]]
    if (ax.type == "predictive") {
        z.axes <- lapply(1:p, function(k, unscaled.X, means, 
            sd, n.int, D, Y, lambda, s, e.vecs, prediction.type) {
            number.points <- num.points
            std.markers <- pretty(unscaled.X[, k], n = n.int[k])
            std.markers <- std.markers[std.markers >= min(unscaled.X[, 
                k]) & std.markers <= max(unscaled.X[, k])]
            interval <- (c(std.markers, min(unscaled.X[, k]), 
                max(unscaled.X[, k])) - means[k])/sd[k]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            axis.vals <- axis.vals[!axis.vals == 0]
            number.points <- length(axis.vals)
            n <- nrow(unscaled.X)
            p <- ncol(unscaled.X)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            mu <- axis.vals
            m <- length(mu)
            ddmu.dnplus1.mat <- t(sapply(unscaled.X[, k], function(x) eval(ddmu.expr(x))))
            A.mat <- t(ddmu.dnplus1.mat) %*% Y[, e.vects] %*% 
                solve(lambda[e.vects, e.vects])
            m.star.vec <- apply(ddmu.dnplus1.mat, 2, sum)/(-n)
            a1.sq.plus.a2.sq <- apply(A.mat, 1, function(a) sum(a^2))
            LL.mat <- cbind(A.mat, m.star.vec, A.mat/sqrt(a1.sq.plus.a2.sq), 
                m.star.vec/sqrt(a1.sq.plus.a2.sq))
            L.a <- function(mu.val, X, k) {
                ddmu.dnplus1.vec <- sapply(X[, k], function(x, 
                  mu) eval(ddmu.expr(x)), mu = mu.val)
                t(ddmu.dnplus1.vec) %*% Y[, e.vects] %*% solve(lambda[e.vects, 
                  e.vects])
            }
            L.m <- function(mu.val, X, k) {
                ddmu.dnplus1.vec <- sapply(X[, k], function(x, 
                  mu) eval(ddmu.expr(x)), mu = mu.val)
                m.star <- sum(ddmu.dnplus1.vec)/(-nrow(X))
                m.star/sum(L.a(mu.val, X, k)^2)
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
                epsilon <- 0.001
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
                    interval.einde)], X = unscaled.X, k)
                  mu.0 <- uit$root
                }
                skep.integrand.1 <- function(mu, X, k) {
                  sapply(mu, function(mu, Xmat, k) {
                    ddmu.dnplus1.vec <- sapply(X[, k], function(x, 
                      mu) eval(ddmu.expr(x)), mu = mu)
                    a.vec <- t(ddmu.dnplus1.vec) %*% Y[, e.vects] %*% 
                      solve(lambda[e.vects, e.vects])
                    wortel <- sqrt(sum(a.vec^2))
                    l1.mu <- a.vec[1]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-nrow(X))
                    d2.dmu2.dnplus1.vec <- sapply(Xmat[, k], 
                      function(x, mu) eval(d2dmu2.expr(x)), mu = mu)
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-nrow(X))
                    ddmu.a1 <- t(Y[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(Y[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l2mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[2] - (ddmu.a2 * 
                      wortel - ddmu.wortel * a.vec[2]) * mster.mu)/(wortel * 
                      a.vec[2]^2)
                    l1.mu * ddmu.mmu.l2mu
                  }, Xmat = X, k = k)
                }
				
                f1.int <- lapply(mu, function(mu, unscaled.X, k) 
				     try(integrate(skep.integrand.1, mu.0, mu, 
					 unscaled.X, k)$value, silent=TRUE), unscaled.X = unscaled.X, k = k)  
				f1.int <- lapply(f1.int, function(x) if(is.null(attr(x,"condition"))) x <- x
                         else x <- attr(x,"condition"))
				f1.int <- unlist(lapply(f1.int, function(x)ifelse(is.numeric(x),x,NA)))
  
                axis.points[, 1] <- LL.mat[, 5] * f1.int
				
                skep.integrand.2 <- function(mu, X, k) {
                  sapply(mu, function(mu, Xmat, k) {
                    ddmu.dnplus1.vec <- sapply(X[, k], function(x, 
                      mu) eval(ddmu.expr(x)), mu = mu)
                    a.vec <- t(ddmu.dnplus1.vec) %*% Y[, e.vects] %*% 
                      solve(lambda[e.vects, e.vects])
                    wortel <- sqrt(sum(a.vec^2))
                    l2.mu <- a.vec[2]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-nrow(X))
                    d2.dmu2.dnplus1.vec <- sapply(Xmat[, k], 
                      function(x, mu) eval(d2dmu2.expr(x)), mu = mu)
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-nrow(X))
                    ddmu.a1 <- t(Y[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(Y[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l1mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[1] - (ddmu.a1 * 
                      wortel - ddmu.wortel * a.vec[1]) * mster.mu)/(wortel * 
                      a.vec[1]^2)
                    l2.mu * ddmu.mmu.l1mu
                  }, Xmat = X, k = k)
                }

				f2.int <- lapply(mu, function(mu, unscaled.X, k) 
				     try(integrate(skep.integrand.2, mu.0, mu, 
					 unscaled.X, k)$value, silent=TRUE), unscaled.X = unscaled.X, k = k)  
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
        }, unscaled.X = unscaled.X, means = means, sd = sd, n.int = n.int, 
            D = D, Y = U, lambda = lambda, s = s, e.vecs = e.vects, 
            prediction.type = prediction.type)
    }
    if (ax.type == "interpolative") {
        z.axes <- lapply(1:p, function(j, unscaled.X, means, 
            sd, n.int, D, Y, lambda.inv, s, dist.func, e.vecs, 
            ax) {
            number.points <- 100
            std.markers <- pretty(unscaled.X[, j], n = n.int[j])
            std.range <- c(min(std.markers), max(std.markers))
            std.markers.min <- std.markers - (std.range[2] - 
                std.range[1])
            std.markers.max <- std.markers + (std.range[2] - 
                std.range[1])
            interval <- (c(std.markers, std.markers.min, std.markers.max) - 
                means[j])/sd[j]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            number.points <- length(axis.vals)
            mat <- matrix(0, ncol = ncol(unscaled.X), nrow = length(axis.vals))
            mat[, j] <- axis.vals
            mat <- rbind(unscaled.X, mat)
            D.nuut <- dist.func(mat)
            dnplus1.mat <- D.nuut[1:n, -(1:n)]
            as.embed <- apply(dnplus1.mat, 2, function(dn1, D, 
                Y, lambda.inv, s) lambda.inv %*% t(Y) %*% (dn1 - 
                D %*% s), D = D, Y = Y, lambda.inv = lambda.inv, 
                s = s)
            axis.points <- cbind(t(as.embed)[, e.vecs], 0, 0)
            axis.points[, 3] <- axis.vals * sd[j] + means[j]
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            if (scale.axes) {
                axis.points[, 1:2] <- scale(axis.points[, 1:2], 
                  center = O.co, scale = FALSE)
                axis.points[, 1:2] <- length(ax) * axis.points[, 
                  1:2]
                axis.points[, 1] <- axis.points[, 1] + O.co[1]
                axis.points[, 2] <- axis.points[, 2] + O.co[2]
            }
            return(axis.points)
        }, unscaled.X = unscaled.X, means = means, sd = sd, n.int = n.int, 
            D = D, Y = U, lambda.inv = lambda.inv, s = s, dist.func = D.mat, 
            e.vecs = e.vects[1:2], ax = ax)
    }
	z.axes <- lapply(z.axes, na.omit)
    uit <- drawbipl.nonlin(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
        p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
        alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
        pch.samples.size = pch.samples.size, specify.bags = NULL, 
        label = label, markers = markers, Title = Title, large.scale = large.scale, 
        specify.classes, Tukey.median = Tukey.median, Z.means.mat = NULL, 
        offset = offset, pos = pos, strepie = line.length, max.num = max.num, 
        c.hull.n = c.hull.n, marker.size = marker.size, label.size = label.size, 
        exp.factor = exp.factor, line.width = line.width, class.vec = class.vec, 
        predictions.sample = predictions.sample, predictions.mean = predictions.mean, 
        ort.lty = ort.lty, straight = straight)
    points(0, 0, pch = 3, cex = 5, col = "light grey")
    if (!is.null(X.new.samples)) {
        lambda.inv.full <- cbind(lambda.inv, 0)
        lambda.inv.full <- rbind(lambda.inv.full, 0)
        ynew.mat <- matrix(NA, nrow = nrow(X.new.samples), ncol = n)
        for (i in 1:nrow(X.new.samples)) {
            dvec <- dvec.fun(X.new.samples[i, ], X)
            ynew.mat[i, ] <- lambda.inv.full %*% t(swd$u %*% 
                diag(swd$d^0.5)) %*% (dvec - apply(D, 1, sum)/n)
        }
        points(x = ynew.mat[, 1], y = ynew.mat[, 2], pch = pch.new, 
            cex = pch.new.size, col = pch.new.cols)
        text(x = ynew.mat[, 1], y = ynew.mat[, 2], labels = pch.new.labels, 
            cex = pch.new.labels.size, col = "black", pos = 1)
    }
    if (!is.null(zoomval)) {
        zoomval <- zoom(zoomval)
        drawbipl.nonlin(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
            alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples.size = pch.samples.size, specify.bags = NULL, 
            label = label, markers = markers, Title = Title, 
            large.scale = large.scale, specify.classes, Tukey.median = Tukey.median, 
            Z.means.mat = NULL, offset = offset, pos = pos, strepie = line.length, 
            max.num = max.num, c.hull.n = c.hull.n, marker.size = marker.size, 
            label.size = label.size, exp.factor = exp.factor, 
            line.width = line.width, class.vec = class.vec, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, ort.lty = ort.lty, 
            straight = straight, zoomval = zoomval)
        points(0, 0, pch = 3, cex = 5, col = "light grey")
    }
    list(D = D, Z = Z, Z.axes = z.axes, e.vals = eigval, quality = fit.quality, 
        predictions = uit$predictions, O.co = O.co, ynew.mat = ynew.mat, 
        B = B)
}
