Genbipl <-
function (X = Remuneration.data.genbipl.2002, G = NULL, X.new.samples = NULL, 
    cont.scale = c("none", "unitVar", "unitSS", "unitRange"), 
    e.vects = 1:2, alpha = 0.95, ax = 1:ncol(X), ax.name.size = 0.65, 
    ax.type = c("predictive", "interpolative"), ax.col = list(ax.col = rep(8, 
        ncol(X)), tickmarker.col = rep(8, ncol(X)), marker.col = rep(1, 
        ncol(X))), c.hull.n = 10, CLPs.plot = FALSE, colours = UBcolours, 
    dist.cont = c("Pythagoras", "Clark", "SqrtL1"), dist.cat = "EMC", 
    exp.factor = 1.2, label = TRUE, label.size = 0.6, line.length = c(1, 
        1), line.type = 1, line.width = 1, markers = TRUE, marker.size = 0.5, 
    max.num = 2500, n.int = rep(5, ncol(X)), offset = rep(0.5, 
        4), ort.lty = 1, parplotmar = rep(3, 4), pch.means = 0:10, 
    pch.means.size = 1, pch.samples = rep(15, 10), pch.samples.size = 1, 
    pos = c("Orthog", "Hor", "Paral"), prediction.regions = FALSE, 
    predictions.sample = NULL, prediction.type = c("normal", 
        "circle"), specify.bags = dimnames(G)[[2]], specify.classes = dimnames(G)[[2]], 
    straight = FALSE, Title = NULL, Tukey.median = FALSE, zoomval = NULL, 
    x.grid = 0.25, y.grid = 0.25, plot.symbol = 20, plot.symbol.size = 1, 
    colours.pred.regions = UBcolours) 
{
    cont.scale <- cont.scale[1]
    if (is.na(match(cont.scale, c("none", "unitVar", "unitSS", 
        "unitRange")))) 
        stop("cont.scale must be one of: 'none','unitVar','unitSS','unitRange' \n")
    dist.cont <- dist.cont[1]
    prediction.type <- prediction.type[1]
    ax.type <- ax.type[1]
    pos <- pos[1]
    old.par <- par(no.readonly = TRUE)
    par(pty = "s", mar = parplotmar)
    if (dist.cont == "Clark") 
        dist.expr <- function(xik, xjk) {
            ((xik - xjk)/(xik + xjk))^2
        }
    if (dist.cont == "Pythagoras") 
        dist.expr <- function(xik, xjk) {
            (xik - xjk)^2
        }
    if (dist.cont == "SqrtL1") 
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
    D.cat.mat <- function(Gklist, numcat, cat.scale) {
        Dk.list <- vector("list", numcat)
        for (k in 1:numcat) {
            nn <- nrow(Gklist[[k]])
            tempmat <- matrix(1, nrow = nn, ncol = nn) - Gklist[[k]] %*% 
                t(Gklist[[k]])
            if (cat.scale == "EMC") 
                tempmat <- tempmat/(sum(tempmat)/(2 * nn))
            Dk.list[[k]] <- -0.5 * (tempmat)
        }
        Dk.list
    }
    dnplus1 <- function(X, x.new) {
        X1 <- X[, !categorical, drop = F]
        X2 <- X[, categorical, drop = F]
        d1 <- apply(X1, 1, function(x, x.new) -0.5 * sum(dist.expr(x, 
            x.new)), x.new = x.new[!categorical, drop = F])
        d2 <- D.cat.mat(rbind(x.new[categorical, drop = F], X2))[1, 
            -1]
        d1 + d2
    }
    ddmu.expr <- function(x) {
        if (dist.cont == "Pythagoras" || dist.cont == "Clark") {
            dist.text <- deparse(dist.expr)[3]
            dist.text <- paste("-0.5", dist.text, sep = "*")
            dist.text <- gsub("xik", x, dist.text)
            dist.text <- gsub("xjk", "mu", dist.text)
            dist.expression <- parse(text = dist.text)
            AA <- D(dist.expression, "mu")
        }
        if (dist.cont == "SqrtL1") 
            AA <- expression(-(sign(x - mu)))
        return(AA)
    }
    d2dmu2.expr <- function(x) {
        first <- ddmu.expr(x)
        D(first, "mu")
    }
    sumlist <- function(listmat) {
        sum <- 0
        for (i in 1:length(listmat)) sum <- sum + listmat[[i]]
        sum
    }
    n <- nrow(X)
    p <- ncol(X)
    categorical <- rep(F, p)
    for (j in 1:p) if (is.factor(X[, j])) 
        categorical[j] <- T
    p2 <- sum(categorical)
    p1 <- p - p2
    unscaled.X <- X
    gems <- rep(0, p1)
    sds <- rep(1, p1)
    if (cont.scale == "unitVar") {
        gems <- apply(X[, !categorical], 2, mean)
        sds <- sqrt(apply(X[, !categorical], 2, var))
        X[, !categorical] <- scale(X[, !categorical])
    }
    if (cont.scale == "unitSS") {
        gems <- apply(X[, !categorical, drop = FALSE], 2, mean)
        sds <- sqrt((n - 1) * apply(X[, !categorical, drop = FALSE], 
            2, var))
        X[, !categorical] <- scale(X[, !categorical, drop = FALSE], 
            center = TRUE, scale = sds)
    }
    if (cont.scale == "unitRange") {
        gems <- apply(X[, !categorical], 2, mean)
        sds <- apply(X[, !categorical], 2, function(x) max(x) - 
            min(x))
        X[, !categorical] <- scale(X[, !categorical], center = TRUE, 
            scale = sds)
    }
    Xmat <- X
    Gmat <- NULL
    Gk.list <- vector("list", p2)
    L <- 0
    for (k in 1:p2) {
        Gk <- indmat(Xmat[, categorical, drop = F][, k])
        Gk.list[[k]] <- Gk
        Gmat <- cbind(Gmat, Gk)
        L <- L + ncol(Gk)
    }
    s <- rep(1/n, n)
    one <- rep(1, n)
    N <- one %*% t(s)
    I <- diag(n)
    Delta1 <- D.mat(Xmat[, !categorical, drop = F])
    Delta2.list <- D.cat.mat(Gk.list, numcat = p2, cat.scale = dist.cat)
    Delta2 <- sumlist(Delta2.list)
    D <- Delta1 + Delta2
    B <- (I - N) %*% D %*% (I - N)
    swd <- svd(B)
    U <- (swd$u %*% diag(swd$d^0.5))[, -n]
    lambda <- t(U) %*% U
    lambda.inv <- diag(ifelse(zapsmall(diag(lambda)) > 0, 1/diag(lambda), 
        0))
    eigval.pos <- sum(zapsmall(swd$d) > 0)
    eigval <- diag(lambda)
    eigval.r <- eigval[e.vects[1:2]]
    fit.predictivity.mat <- NULL
    fit.predictivity <- NULL
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
        "%")
    fit.adequacy <- NULL
    CLPs <- NULL
    Z <- U[, e.vects[1:2]]
    if (is.null(G)) 
        J <- 0
    else J <- ncol(G)
    Z <- data.frame(Z, pch.sampl = pch.samples[1], colr = as.character(colours[1]), 
        line.type = line.type[1], stringsAsFactors = FALSE)
    if (J > 0) {
        for (j in 1:J) {
            Z[G[, j] == 1, 4] <- as.character(colours[j])
            Z[G[, j] == 1, 3] <- pch.samples[j]
            Z[G[, j] == 1, 5] <- line.type[j]
        }
    }
    rownames(Z) <- rownames(Xmat)
    if (ax.type == "predictive") {
        z.axes <- lapply(1:p1, function(k, Xmat, unscaled.X, 
            means, sd, n.int, D, G, lambda, s, e.vecs, prediction.type) {
            number.points <- 100
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
            ddmu.dnplus1.mat <- t(sapply(Xmat[, k], function(x) eval(ddmu.expr(x))))
            A.mat <- t(ddmu.dnplus1.mat) %*% G[, e.vects] %*% 
                solve(lambda[e.vects, e.vects])
            m.star.vec <- apply(ddmu.dnplus1.mat, 2, sum)/(-n)
            a1.sq.plus.a2.sq <- apply(A.mat, 1, function(a) sum(a^2))
            LL.mat <- cbind(A.mat, m.star.vec, A.mat/sqrt(a1.sq.plus.a2.sq), 
                m.star.vec/sqrt(a1.sq.plus.a2.sq))
            L.a <- function(mu.val, X, k) {
                ddmu.dnplus1.vec <- sapply(X[, k], function(x, 
                  mu) eval(ddmu.expr(x)), mu = mu.val)
                t(ddmu.dnplus1.vec) %*% G[, e.vects] %*% solve(lambda[e.vects, 
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
                }, unscaled.X = Xmat, k = k)
                xsi.k <- solve(lambda[e.vects, e.vects]) %*% 
                  t(G[, e.vects]) %*% (mat - apply(D, 1, sum)/n)
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
                    interval.einde)], X = Xmat, k)
                  mu.0 <- uit$root
                }
                skep.integrand.1 <- function(mu, XX, k) {
                  sapply(mu, function(mu, Xmat, k) {
                    ddmu.dnplus1.vec <- sapply(XX[, k], function(x, 
                      mu) eval(ddmu.expr(x)), mu = mu)
                    a.vec <- t(ddmu.dnplus1.vec) %*% G[, e.vects] %*% 
                      solve(lambda[e.vects, e.vects])
                    wortel <- sqrt(sum(a.vec^2))
                    l1.mu <- a.vec[1]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-nrow(XX))
                    d2.dmu2.dnplus1.vec <- sapply(Xmat[, k], 
                      function(x, mu) eval(d2dmu2.expr(x)), mu = mu)
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-nrow(XX))
                    ddmu.a1 <- t(G[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(G[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l2mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[2] - (ddmu.a2 * 
                      wortel - ddmu.wortel * a.vec[2]) * mster.mu)/(wortel * 
                      a.vec[2]^2)
                    l1.mu * ddmu.mmu.l2mu
                  }, Xmat = Xmat, k = k)
                }
                f1.int <- sapply(mu, function(mu, Xmat, k) my.integrate(skep.integrand.1, 
                  mu.0, mu, Xmat, k)$value, Xmat = Xmat, k = k)
                axis.points[, 1] <- LL.mat[, 5] * f1.int
                skep.integrand.2 <- function(mu, XX, k) {
                  sapply(mu, function(mu, Xmat, k) {
                    ddmu.dnplus1.vec <- sapply(XX[, k], function(x, 
                      mu) eval(ddmu.expr(x)), mu = mu)
                    a.vec <- t(ddmu.dnplus1.vec) %*% G[, e.vects] %*% 
                      solve(lambda[e.vects, e.vects])
                    wortel <- sqrt(sum(a.vec^2))
                    l2.mu <- a.vec[2]/wortel
                    mster.mu <- sum(ddmu.dnplus1.vec)/(-nrow(XX))
                    d2.dmu2.dnplus1.vec <- sapply(Xmat[, k], 
                      function(x, mu) eval(d2dmu2.expr(x)), mu = mu)
                    ddmu.mster.mu <- sum(d2.dmu2.dnplus1.vec)/(-nrow(XX))
                    ddmu.a1 <- t(G[, 1]) %*% d2.dmu2.dnplus1.vec/lambda[1, 
                      1]
                    ddmu.a2 <- t(G[, 2]) %*% d2.dmu2.dnplus1.vec/lambda[2, 
                      2]
                    ddmu.wortel <- (a.vec[1] * ddmu.a1 + a.vec[2] * 
                      ddmu.a2)/wortel
                    ddmu.mmu.l1mu <- ((ddmu.mster.mu * wortel - 
                      ddmu.wortel * mster.mu) * a.vec[1] - (ddmu.a1 * 
                      wortel - ddmu.wortel * a.vec[1]) * mster.mu)/(wortel * 
                      a.vec[1]^2)
                    l2.mu * ddmu.mmu.l1mu
                  }, Xmat = Xmat, k = k)
                }
                f2.int <- sapply(mu, function(mu, Xmat, k) my.integrate(skep.integrand.2, 
                  mu.0, mu, Xmat, k)$value, Xmat = Xmat, k = k)
                axis.points[, 2] <- LL.mat[, 4] * f2.int
            }
            axis.points[, 3] <- zapsmall(axis.vals * sd[k] + 
                means[k])
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, Xmat = Xmat[, !categorical, drop = F], unscaled.X = unscaled.X[, 
            !categorical, drop = F], means = gems, sd = sds, 
            n.int = n.int, D = D, G = U, lambda = lambda, s = s, 
            e.vecs = e.vects, prediction.type = prediction.type)
        CLPs <- lapply(1:p2, function(k, Xmat, D, G, lambda.inv, 
            s, e.vecs) {
            axis.vals <- levels(X[, categorical, drop = F][, 
                k])
            dnplus1.mat <- sapply(axis.vals, function(mu, Xmat) {
                X.pseudo <- Xmat
                X.pseudo[, k] <- mu
                X2.both <- rbind(Xmat, X.pseudo)
                Gk.list <- vector("list", p2)
                for (j in 1:p2) {
                  Gk.list[[j]] <- indmat(X2.both[, j])
                }
                D2.both.list <- D.cat.mat(Gk.list, numcat = p2, 
                  cat.scale = dist.cat)
                D2.both <- sumlist(D2.both.list)
                D11.1 <- Delta1
                D22.2 <- D2.both[-(1:n), -(1:n)]
                D12.2 <- D2.both[1:n, -(1:n)]
                sum(D11.1 + D22.2) * matrix(1, nrow = n, ncol = 1)/(-2 * 
                  n) + (D11.1 + D12.2) %*% matrix(1, nrow = n, 
                  ncol = 1)/n
            }, Xmat = Xmat)
            CLPs <- t(apply(dnplus1.mat, 2, function(dn1, D, 
                G, lambda.inv, s) (lambda.inv[1:e.vecs, 1:e.vecs]) %*% 
                t(G[, 1:e.vecs]) %*% (dn1 - D %*% s), D = D, 
                G = G, lambda.inv = lambda.inv, s = s))
            return(CLPs)
        }, Xmat = Xmat[, categorical, drop = F], D = D, G = U, 
            lambda.inv = lambda.inv, s = s, e.vecs = eigval.pos)
    }
    if (ax.type == "interpolative") {
        z.axes <- lapply(1:p1, function(k, Xmat, unscaled.X, 
            means, sd, n.int, D, G, lambda.inv, s, dist.func, 
            e.vecs, ax) {
            number.points <- 100
            std.markers <- pretty(unscaled.X[, k], n = n.int[k])
            std.markers <- std.markers[std.markers >= min(unscaled.X[, 
                k]) & std.markers <= max(unscaled.X[, k])]
            interval <- (c(std.markers, min(unscaled.X[, k]), 
                max(unscaled.X[, k])) - means[k])/sd[k]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            axis.vals <- zapsmall(axis.vals)
            if (dist.cont == "Clark") 
                axis.vals <- axis.vals[!axis.vals == 0]
            number.points <- length(axis.vals)
            dnplus1.mat <- sapply(axis.vals, function(mu, Xmat, 
                means, sd, k) {
                X.pseudo <- Xmat
                X.pseudo[, k] <- mu
                X1.both <- rbind(Xmat, X.pseudo)
                D1.both <- D.mat(X1.both)
                D22.1 <- D1.both[-(1:n), -(1:n)]
                D11.2 <- Delta2
                D12.1 <- D1.both[1:n, -(1:n)]
                sum(D22.1 + D11.2) * matrix(1, nrow = n, ncol = 1)/(-2 * 
                  n) + (D12.1 + D11.2) %*% matrix(1, nrow = n, 
                  ncol = 1)/n
            }, Xmat = Xmat, means = means, sd = sds, k = k)
            as.embed <- apply(dnplus1.mat, 2, function(dn1, D, 
                G, lambda.inv, s) lambda.inv %*% t(G) %*% (dn1 - 
                D %*% s), D = D, G = G, lambda.inv = lambda.inv, 
                s = s)
            axis.points <- cbind(t(as.embed)[, e.vecs], 0, 0)
            axis.points[, 3] <- axis.vals * sd[k] + means[k]
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, Xmat = Xmat[, !categorical, drop = F], unscaled.X = unscaled.X[, 
            !categorical, drop = F], means = gems, sd = sds, 
            n.int = n.int, D = D, G = U, lambda.inv = lambda.inv, 
            s = s, dist.func = D.mat, e.vecs = e.vects[1:2], 
            ax = ax)
        CLPs <- lapply(1:p2, function(k, Xmat, D, G, lambda.inv, 
            s, e.vecs) {
            axis.vals <- levels(X.mat[, categorical, drop = F][, 
                k])
            dnplus1.mat <- sapply(axis.vals, function(mu, Xmat) {
                X.pseudo <- Xmat
                X.pseudo[, k] <- mu
                X2.both <- rbind(Xmat, X.pseudo)
                Gk.list <- vector("list", p2)
                for (j in 1:p2) {
                  Gk.list[[j]] <- indmat(X2.both[, j])
                }
                D2.both.list <- D.cat.mat(Gk.list, numcat = p2, 
                  cat.scale = dist.cat)
                D2.both <- sumlist(D2.both.list)
                D11.1 <- Delta1
                D22.2 <- D2.both[-(1:n), -(1:n)]
                D12.2 <- D2.both[1:n, -(1:n)]
                sum(D11.1 + D22.2) * matrix(1, nrow = n, ncol = 1)/(-2 * 
                  n) + (D11.1 + D12.2) %*% matrix(1, nrow = n, 
                  ncol = 1)/n
            }, Xmat = Xmat)
            CLPs <- t(apply(dnplus1.mat, 2, function(dn1, D, 
                G, lambda.inv, s) (lambda.inv[1:e.vecs, 1:e.vecs]) %*% 
                t(G[, 1:e.vecs]) %*% (dn1 - D %*% s), D = D, 
                G = G, lambda.inv = lambda.inv, s = s))
            return(CLPs)
        }, Xmat = Xmat[, categorical, drop = F], D = D, G = U, 
            lambda.inv = lambda.inv, s = s, e.vecs = eigval.pos)
    }
    uit <- drawbipl.genbipl(Z, z.axes, z.axes.names = dimnames(Xmat)[[2]], 
        p = p1, ax = ax[1:p1], ax.col = ax.col, ax.name.size = ax.name.size, 
        alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
        pch.samples.size = pch.samples.size, colours = colours, 
        CLPs.plot = CLPs.plot, specify.bags = NULL, label = label, 
        markers = markers, Title = Title, specify.classes = specify.classes, 
        Tukey.median = Tukey.median, Z.means.mat = NULL, offset = offset, 
        pos = pos, strepie = line.length, max.num = max.num, 
        c.hull.n = c.hull.n, marker.size = marker.size, label.size = label.size, 
        exp.factor = exp.factor, line.width = line.width, class.vec = class.vec, 
        predictions.sample = predictions.sample, predictions.mean = NULL, 
        ort.lty = ort.lty, parplotmar = parplotmar, straight = straight, 
        CLPs = CLPs, prediction.regions = prediction.regions, 
        x.grid = x.grid, y.grid = y.grid, plot.symbol = plot.symbol, 
        plot.symbol.size = plot.symbol.size, colours.pred.regions = colours.pred.regions)
    points(0, 0, pch = 3, cex = 5, col = "light grey")
    if (!is.null(zoomval)) {
        zoomval <- zoom(zoomval)
        s.drawbipl(Z, z.axes, z.axes.names = dimnames(Xmat)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
            alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples.size = pch.samples.size, colours = colours, 
            specify.bags = NULL, label = label, markers = markers, 
            Title = Title, large.scale = large.scale, specify.classes, 
            Tukey.median = Tukey.median, Z.means.mat = NULL, 
            offset = offset, pos = pos, strepie = line.length, 
            max.num = max.num, c.hull.n = c.hull.n, marker.size = marker.size, 
            label.size = label.size, exp.factor = exp.factor, 
            line.width = line.width, class.vec = class.vec, predictions.sample = predictions.sample, 
            predictions.mean = NULL, ort.lty = ort.lty, straight = straight, 
            zoomval = zoomval)
        points(0, 0, pch = 3, cex = 5, col = "light grey")
    }
    list(Z = Z, Z.axes = z.axes, e.vals = eigval, D.cont = Delta1, 
        D.cat = Delta2, normalized.X.cont = Xmat[, !categorical], 
        centred.vec = gems, scaling.vec = sds, CLPs = CLPs, usr = uit$usr, 
        predictions = uit$predictions, category.predictions.list = uit$category.predictions.list)
}
