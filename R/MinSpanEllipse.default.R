MinSpanEllipse.default <-
function (x, asp, clus, diss = FALSE, cor = TRUE, stand = FALSE, 
    lines = 2, shade = FALSE, color = FALSE, labels = 0, plotchar = TRUE, 
    col.p = "dark green", col.txt = col.p, col.clus = if (color) c(2, 
        4, 6, 3) else 5, span = TRUE, xlim = NULL, ylim = NULL, 
    main = "", sub = "", verbose = getOption("verbose"), pch.points = 16, 
    ...) 
{
    names <- names(x)
    force(main)
    if (is.data.frame(x)) 
        x <- data.matrix(x)
    if (!is.numeric(x)) 
        stop("x is not numeric")
    if (diss) {
        if (any(is.na(x))) 
            stop("NA-values are not allowed in dist-like `x'.")
        if (inherits(x, "dist")) {
            n <- attr(x, "Size")
            labels1 <- attr(x, "Labels")
        }
        else {
            siz <- sizeDiss(x)
            if (is.na(siz)) {
                if ((n <- nrow(x)) != ncol(x)) 
                  stop("Distances must be result of dist or a square matrix.")
                if (all.equal(x, t(x)) != TRUE) 
                  stop("the square matrix is not symmetric.")
                labels1 <- dimnames(x)[[1]]
            }
            else {
                if (!is.vector(x)) {
                  labels1 <- attr(x, "Labels")
                  x <- as.matrix(x)
                  if ((n <- nrow(x)) == ncol(x) && all.equal(x, 
                    t(x)) == TRUE) {
                    labels1 <- dimnames(x)[[1]]
                  }
                  else {
                    warning(">>>>> funny case in clusplot.default() -- please report!\n")
                    if (is.null(labels1)) 
                      labels1 <- 1:sizeDiss(x)
                    attr(x, "Size") <- sizeDiss(x)
                  }
                }
                else {
                  attr(x, "Size") <- n <- siz
                  labels1 <- 1:n
                }
            }
        }
        if (is.null(labels1)) 
            labels1 <- 1:n
        if (paste(R.version$major, R.version$minor, sep = ".") < 
            1.5) {
            cmdscale <- function(d, k = 2, add = TRUE, ...) {
                if (any(is.na(d))) 
                  stop("NA values not allowed in d")
                if (is.null(n <- attr(d, "Size"))) {
                  d <- as.matrix(d)
                  x <- d^2
                  if ((n <- nrow(x)) != ncol(x)) 
                    stop("Distances must be result of dist or a square matrix")
                }
                else {
                  x <- matrix(0, n, n)
                  if (add) 
                    d0 <- x
                  x[row(x) > col(x)] <- d^2
                  x <- x + t(x)
                  if (add) {
                    d0[row(x) > col(x)] <- d
                    d <- d0 + t(d0)
                  }
                }
                storage.mode(x) <- "double"
                x <- .C("dblcen", x = x, as.integer(n), PACKAGE = "mva")$x
                if (add) {
                  i2 <- n + (i <- 1:n)
                  Z <- matrix(0, 2 * n, 2 * n)
                  Z[cbind(i2, i)] <- -1
                  Z[i, i2] <- -x
                  Z[i2, i2] <- .C("dblcen", x = 2 * d, as.integer(n), 
                    PACKAGE = "mva")$x
                  e <- eigen(Z, symmetric = FALSE, only.val = TRUE)$values
                  add.c <- max(Re(e))
                  x <- matrix(double(n * n), n, n)
                  non.diag <- row(d) != col(d)
                  x[non.diag] <- (d[non.diag] + add.c)^2
                }
                e <- eigen(-x/2, symmetric = TRUE)
                ev <- e$values[1:k]
                points <- e$vectors[, 1:k] %*% diag(sqrt(ev), 
                  k)
                rn <- if (is.matrix(d)) 
                  rownames(d)
                else names(d)
                dimnames(points) <- list(rn, NULL)
                evalus <- e$values[-n]
                list(points = points, eig = ev, ac = if (add) add.c else 0, 
                  GOF = sum(ev)/c(sum(abs(evalus)), sum(evalus[evalus > 
                    0])))
            }
        }
        x1 <- cmdscale(x, k = 2, eig = TRUE, add = TRUE)
        if (x1$ac < 0) 
            x1 <- cmdscale(x, k = 2, eig = TRUE)
        var.dec <- x1$GOF[2]
        x1 <- x1$points
    }
    else {
        if (!is.matrix(x)) 
            stop("x is not a data matrix")
        if (any(is.na(x))) {
            y <- is.na(x)
            if (any(apply(y, 1, all))) 
                stop("one or more objects contain only missing values")
            if (any(apply(y, 2, all))) 
                stop("one or more variables contain only missing values")
            x <- apply(x, 2, function(x) {
                x[is.na(x)] <- median(x, na.rm = TRUE)
                x
            })
            cat("Missing values were displaced by the median of the corresponding variable(s)\n")
        }
        n <- nrow(x)
        labels1 <- dimnames(x)[[1]]
        if (is.null(labels1)) 
            labels1 <- 1:n
        if (ncol(x) == 1) {
            hulp <- rep(0, length(x))
            x1 <- matrix(c(t(x), hulp), ncol = 2)
            var.dec <- 1
        }
        else {
            prim.pr <- princomp(x, scores = TRUE, cor = ncol(x) != 
                2)
            var.dec <- cumsum(prim.pr$sdev^2/sum(prim.pr$sdev^2))[2]
            x1 <- prim.pr$scores
            x1 <- cbind(x1[, 1], x1[, 2])
        }
    }
    x1 <- x
    clus <- as.vector(clus)
    if (length(clus) != n) 
        stop("The clustering vector has not the good length")
    clus <- as.factor(clus)
    if (any(is.na(clus))) 
        stop("NA-values are not allowed in clustering vector")
    if (stand) 
        x1 <- scale(x1)
    rangx <- range(x1[, 1])
    rangy <- range(x1[, 2])
    minx <- rangx[1]
    maxx <- rangx[2]
    miny <- rangy[1]
    maxy <- rangy[2]
    levclus <- levels(clus)
    nC <- length(levclus)
    z <- A <- vector("list", nC)
    maxima <- loc <- matrix(0, nrow = nC, ncol = 2)
    d2 <- verhoud <- numeric(nC)
    verhouding <- 0
    num3 <- 90
    num6 <- 70
    for (i in 1:nC) {
        x <- x1[clus == levclus[i], , drop = FALSE]
        aantal <- nrow(x)
        cov <- var(if (aantal == 1) {
            if (verbose) 
                cat("cluster", i, " has only one observation ..\n")
            rbind(x, c(0, 0))
        }
        else x)
        x.1 <- range(x[, 1])
        y.1 <- range(x[, 2])
        notrank2 <- qr(cov, tol = 0.001)$rank != 2
        if (!span && notrank2) {
            d2[i] <- 1
            if ((abs(diff(x.1)) > (diff(rangx)/70)) || (abs(diff(y.1)) > 
                (diff(rangy)/50))) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + 
                  diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 + (loc[i, 2] - 
                  y.1[1])^2)
                a <- a + 0.05 * a
                num2 <- 40
                if (abs(diff(x.1)) > diff(rangx)/70) {
                  ind1 <- which.max(x[, 1])
                  ind2 <- which.min(x[, 1])
                  q <- atan((x[ind1, 2] - x[ind2, 2])/(x[ind1, 
                    1] - x[ind2, 1]))
                  b <- if (diff(rangy) == 0) 
                    1
                  else if (abs(diff(y.1)) > diff(rangy)/50) 
                    diff(y.1)/10
                  else diff(rangy)/num2
                }
                else {
                  b <- if (diff(rangx) == 0) 
                    1
                  else diff(rangx)/num2
                  q <- pi/2
                }
                D <- diag(c(a^2, b^2))
                R <- rbind(c(cos(q), -sin(q)), c(sin(q), cos(q)))
                A[[i]] <- (R %*% D) %*% t(R)
            }
            else {
                a <- diff(rangx)/num3
                b <- diff(rangy)/num6
                if (a == 0) 
                  a <- 1
                if (b == 0) 
                  b <- 1
                A[[i]] <- diag(c(a^2, b^2))
                loc[i, ] <- x[1, ]
            }
            oppervlak <- pi * a * b
        }
        else if (span && notrank2) {
            d2[i] <- 1
            if (sum(x[, 1] != x[1, 1]) != 0 || sum(x[, 2] != 
                x[1, 2]) != 0) {
                loc[i, ] <- c(x.1[1] + diff(x.1)/2, y.1[1] + 
                  diff(y.1)/2)
                a <- sqrt((loc[i, 1] - x.1[1])^2 + (loc[i, 2] - 
                  y.1[1])^2)
                if (any(x[, 1] != x[1, 1])) {
                  ind1 <- which.max(x[, 1])
                  ind2 <- which.min(x[, 1])
                  q <- atan((x[ind1, 2] - x[ind2, 2])/(x[ind1, 
                    1] - x[ind2, 1]))
                }
                else {
                  q <- pi/2
                }
                b <- 1e-07
                D <- diag(c(a^2, b^2))
                R <- rbind(c(cos(q), -sin(q)), c(sin(q), cos(q)))
                A[[i]] <- (R %*% D) %*% t(R)
            }
            else {
                a <- diff(rangx)/num3
                b <- diff(rangy)/num6
                if (a == 0) 
                  a <- 1
                if (b == 0) 
                  b <- 1
                A[[i]] <- diag(c(a^2, b^2))
                loc[i, ] <- x[1, ]
            }
            oppervlak <- pi * a * b
        }
        else {
            if (!span) {
                loc[i, ] <- apply(x, 2, mean)
                d2[i] <- max(mahalanobis(x, loc[i, ], cov))
            }
            else {
                if (verbose) 
                  cat("span & rank2 : calling \"spannel\" ..\n")
                k <- as.integer(2)
                res <- .C("spannel", aantal, ndep = k, dat = cbind(1, 
                  x), sqdist = double(aantal), l1 = double((k + 
                  1)^2), double(k), double(k), prob = double(aantal), 
                  double(k + 1), eps = (0.01), maxit = as.integer(5000), 
                  ierr = integer(1), PACKAGE = "cluster")
                if (res$ierr != 0) 
                  cat("Error in Fortran routine for the spanning ellipsoid,", 
                    "\n rank problem??\n", sep = "")
                cov <- cov.wt(x, res$prob)
                loc[i, ] <- cov$center
                cov <- cov$cov * (1 - sum(cov$wt^2))
                d2[i] <- weighted.mean(res$sqdist, res$prob)
                if (verbose) 
                  cat("ellipse( A= (", format(cov[1, ]), "*", 
                    format(cov[2, 2]), "),\n\td2=", format(d2[i]), 
                    ", loc[]=", format(loc[i, ]), ")\n")
            }
            A[[i]] <- cov
            oppervlak <- pi * d2[i] * sqrt(cov[1, 1] * cov[2, 
                2] - cov[1, 2]^2)
        }
        z[[i]] <- ellipsoidPoints(A[[i]], d2[i], loc[i, ], n = 201)
        maxima[i, ] <- z[[i]][201, ]
        rx <- range(z[[i]][, 1])
        ry <- range(z[[i]][, 2])
        minx <- min(minx, rx[1])
        maxx <- max(maxx, rx[2])
        miny <- min(miny, ry[1])
        maxy <- max(maxy, ry[2])
        verhoud[i] <- aantal/oppervlak
        if (verhoud[i] < 1e+07) 
            verhouding <- verhouding + verhoud[i]
    }
    if (verhouding == 0) 
        verhouding <- 1
    density <- 3 + (verhoud * 37)/verhouding
    density[density > 41] <- 41
    if (span) {
        if (rangx[1] == rangx[2]) {
            minx <- x1[1, 1] - 1
            maxx <- x1[1, 1] + 1
        }
        if (rangy[1] == rangy[2]) {
            miny <- x1[1, 2] - 1
            maxy <- x1[1, 2] + 1
        }
    }
    if (is.null(xlim)) 
        xlim <- c(minx, maxx)
    if (is.null(ylim)) 
        ylim <- c(miny, maxy)
    if (length(col.p) < n) 
        col.p <- rep(col.p, length = n)
    plot(x1[, 1], x1[, 2], asp = asp, xlim = xlim, ylim = ylim, 
        xlab = names[1], ylab = names[2], main = main, type = if (plotchar) 
            "n"
        else "p", col = col.p)
    if (!is.null(sub) && !is.na(sub) && nchar(sub) > 0) 
        title(sub = sub, adj = 0)
    if (color) {
        if (length(col.clus) < min(4, nC)) 
            stop("`col.clus' should have length 4 when color is TRUE")
        i.verh <- order(verhoud)
        jInd <- if (nC > 4) 
            pam(verhoud[i.verh], 4)$clustering
        else 1:nC
        for (i in 1:nC) {
            k <- i.verh[i]
            polygon(z[[k]], density = if (shade) 
                density[k]
            else 0, col = col.clus[jInd[i]], ...)
        }
        col.clus <- col.clus[jInd][order(i.verh)]
    }
    else {
        for (i in 1:nC) polygon(z[[i]], density = if (shade) 
            density[i]
        else 0, col = col.clus, ...)
    }
    if (plotchar) {
        karakter <- pch.points
        for (i in 1:nC) {
            iC <- clus == levclus[i]
            x <- x1[iC, , drop = FALSE]
            il <- 1 + (i - 1)%%19
            points(x[, 1], x[, 2], pch = pch.points, col = col.p[iC], 
                ...)
        }
    }
    if ((lines == 1 || lines == 2) && nC > 1) {
        clas.snijpunt <- function(x, loc, m, n, p) {
            if (loc[n, m] <= x[1, m] && x[1, m] <= loc[p, m]) 
                x[1, ]
            else if (loc[n, m] <= x[2, m] && x[2, m] <= loc[p, 
                m]) 
                x[2, ]
            else NA
        }
        coord.snijp1 <- function(x, gemid) x[2, 2] - 2 * x[1, 
            2] * gemid + x[1, 1] * gemid^2
        coord.snijp2 <- function(x, d2, y) ((x[1, 1] * x[2, 2] - 
            x[1, 2]^2) * d2)/y
        coord.snijp3 <- function(xx, y, gemid) {
            sy <- sqrt(y)
            sy <- c(sy, -sy)
            cbind(xx[1] + sy, xx[2] + gemid * sy)
        }
        afstand <- matrix(0, ncol = nC, nrow = nC)
        for (i in 1:(nC - 1)) {
            for (j in (i + 1):nC) {
                gemid <- (loc[j, 2] - loc[i, 2])/(loc[j, 1] - 
                  loc[i, 1])
                s0 <- coord.snijp1(A[[i]], gemid)
                b0 <- coord.snijp2(A[[i]], d2[i], s0)
                snijp.1 <- coord.snijp3(loc[i, ], y = b0, gemid)
                s1 <- coord.snijp1(A[[j]], gemid)
                b1 <- coord.snijp2(A[[j]], d2[j], s1)
                snijp.2 <- coord.snijp3(loc[j, ], y = b1, gemid)
                if (loc[i, 1] != loc[j, 1]) {
                  if (loc[i, 1] < loc[j, 1]) {
                    punt.1 <- clas.snijpunt(snijp.1, loc, 1, 
                      i, j)
                    punt.2 <- clas.snijpunt(snijp.2, loc, 1, 
                      i, j)
                  }
                  else {
                    punt.1 <- clas.snijpunt(snijp.1, loc, 1, 
                      j, i)
                    punt.2 <- clas.snijpunt(snijp.2, loc, 1, 
                      j, i)
                  }
                }
                else {
                  if (loc[i, 2] < loc[j, 2]) {
                    punt.1 <- clas.snijpunt(snijp.1, loc, 2, 
                      i, j)
                    punt.2 <- clas.snijpunt(snijp.2, loc, 2, 
                      i, j)
                  }
                  else {
                    punt.1 <- clas.snijpunt(snijp.1, loc, 2, 
                      j, i)
                    punt.2 <- clas.snijpunt(snijp.2, loc, 2, 
                      j, i)
                  }
                }
                if (is.na(punt.1[1]) || is.na(punt.2[1]) || (sqrt((punt.1[1] - 
                  loc[i, 1])^2 + (punt.1[2] - loc[i, 2])^2) + 
                  sqrt((punt.2[1] - loc[j, 1])^2 + (punt.2[2] - 
                    loc[j, 2])^2)) > sqrt((loc[j, 1] - loc[i, 
                  1])^2 + (loc[j, 2] - loc[i, 2])^2)) {
                  afstand[i, j] <- NA
                }
                else if (lines == 1) {
                  afstand[i, j] <- sqrt((loc[i, 1] - loc[j, 1])^2 + 
                    (loc[i, 2] - loc[j, 2])^2)
                  segments(loc[i, 1], loc[i, 2], loc[j, 1], loc[j, 
                    2], col = 6, ...)
                }
                else {
                  afstand[i, j] <- sqrt((punt.1[1] - punt.2[1])^2 + 
                    (punt.1[2] - punt.2[2])^2)
                  segments(punt.1[1], punt.1[2], punt.2[1], punt.2[2], 
                    col = 6, ...)
                }
            }
        }
        afstand <- t(afstand) + afstand
    }
    else afstand <- NULL
    if (labels) {
        if (labels == 1) {
            for (i in 1:nC) {
                m <- nrow(z[[i]])
                ni <- length(ii <- seq(1, m, by = max(1, m%/%40)))
                x1 <- rbind(x1, z[[i]][ii, ])
                labels1 <- c(labels1, rep(levclus[i], ni))
            }
            identify(x1, labels = labels1, col = col.txt[1])
        }
        else {
            Stext <- function(xy, labs, ...) {
                xy[, 1] <- xy[, 1] + (maxx - minx)/130
                xy[, 2] <- xy[, 2] + (maxy - miny)/50
                text(xy, labels = labs, ...)
            }
            if (labels == 3 || labels == 2) 
                Stext(x1, labels1, col = col.txt)
            if (labels %in% c(2, 4, 5)) 
                Stext(maxima, levclus, font = 4, col = col.clus)
            if (labels == 5) 
                identify(x1, labels = labels1, col = col.txt[1])
        }
    }
    density[density == 41] <- NA
    invisible(list(Distances = afstand, Shading = density))
}
