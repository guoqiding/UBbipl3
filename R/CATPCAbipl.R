CATPCAbipl <-
function (Xcat = MCA.Table.1.data[, -1], factor.type = rep("nom", 
    ncol(Xcat)), ax = 1:ncol(Xcat), ax.name.size = 0.8, alpha = 95, 
    boxtype = c("n", "o", "l", "7", "c", "u", "]"), drawbagplots = NULL, 
    calibration.pch = 16, calibration.col = "black", calibration.size = 1, 
    calibration.label.size = 1, calibration.label.col = "red", 
    calibration.label.pos = rep(1, ncol(Xcat)), calibration.label.offset = rep(0, 
        ncol(Xcat)), class.vec = Remuneration.cat.data.2002[, 
        3], class.cols = rep(2, length(levels(class.vec))), class.pch = rep(16, 
        length(levels(class.vec))), c.hull.n = 10, epsilon = 1e-06, 
    exp.factor = 1.2, ID = 1:nrow(Xcat), label.samples = FALSE, 
    max.num = 2500, line.type.bags = rep(2, length(levels(class.vec))), 
    line.width = 2, nom.col = 3:12, r = 2, offset = rep(0.5, 
        4), ord.col = rep(5, 10), orthog.transx = rep(0, ncol(Xcat)), 
    orthog.transy = rep(0, ncol(Xcat)), ort.lty = 1, parplotmar = rep(3, 
        4), plot.samples = NULL, pos = c("Orthog", "Hor", "Paral"), 
    predict.sample = NULL, reverse = FALSE, samples.col = "green", 
    samples.pch = 16, samples.size = 0.75, sample.label.pos = 1, 
    sample.label.col = "black", sample.label.offset = 0, sample.label.size = 0.75, 
    select.origin = FALSE, specify.bags = levels(Remuneration.cat.data.2002[, 
        3]), w.factor = 2, z.score.graph = c(3, 2)) 
{
    boxtype <- boxtype[1]
    .bags.plot <- function(Z, specify.bags, class.vec, class.pch, 
        class.cols, usr, alpha1, line.type.bags, max.num, line.width, 
        Tukey.median = FALSE, c.hull.n) {
        classes <- levels(class.vec)[match(specify.bags, levels(class.vec))]
        legend.labs.bags <- classes
        for (j in classes) {
            Z.class <- Z[class.vec == j, ]
            chr <- class.pch[match(j, levels(class.vec))]
            colr <- class.cols[match(j, levels(class.vec))]
            lty <- levels(line.type.bags)[match(j, levels(class.vec))]
            flush.console()
            cat(paste("bag for class ", j, " with ", nrow(Z.class), 
                " samples in class ", j, sep = ""), "\n")
            if (nrow(Z.class) > max.num) 
                Z.class <- Z.class[sample(1:nrow(Z.class), max.num), 
                  ]
            x <- Z.class[, 1]
            y <- Z.class[, 2]
            if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
                if (!is.numeric(x)) 
                  stop(message = "x is not a numeric dataframe or vector.")
            }
            if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
                if ((!is.data.frame(x) && !is.numeric(x)) || 
                  (!all(sapply(x, data.class) == "numeric"))) 
                  stop(message = "x is not a numeric dataframe or vector.")
            }
            x <- as.matrix(x)
            if (dim(x)[2] != 1) 
                stop(message = "x is not a vector.")
            if (is.vector(y) || (is.matrix(y) && !is.data.frame(y))) {
                if (!is.numeric(y)) 
                  stop(message = "y is not a numeric dataframe or vector.")
            }
            if ((!is.matrix(y) && !is.vector(y)) || is.data.frame(y)) {
                if ((!is.data.frame(y) && !is.numeric(y)) || 
                  (!all(sapply(y, data.class) == "numeric"))) 
                  stop(message = "y is not a numeric dataframe or vector.")
            }
            y <- as.matrix(y)
            if (dim(y)[2] != 1) 
                stop(message = "y is not a vector.")
            if (nrow(x) != nrow(y)) 
                stop(message = "x and y should have the same length!")
            na.x <- !is.finite(x)
            na.y <- !is.finite(y)
            ok <- !(na.x | na.y)
            x <- x[ok, , drop = FALSE]
            y <- y[ok, , drop = FALSE]
            n <- nrow(x)
            if (length(x) == 0) 
                stop(message = "All observations have missing values")
            if (n == 1) 
                stop(message = "The sample size should be at least two!")
            dimny <- dimnames(y)[[1]]
            if (length(dimny) == 0) 
                dimny <- 1:n
            if (n < c.hull.n) {
                hull <- chull(x, y)
                points(x, y, pch = as.vector(chr)[1], col = as.vector(colr)[1])
                polygon(x[hull], y[hull], density = 0, col = as.vector(colr)[1], 
                  lwd = line.width, lty = as.vector(lty)[1])
                count <- legend.labs.bags == j
                legend.labs.bags[count] <- paste(legend.labs.bags[count], 
                  ". (C hull)", sep = "")
                warning(paste("Samples too few for alpha bag for class", 
                  j, ". Convex hull constructed with all sample points shown. "))
            }
            else {
                storage.mode(x) <- "double"
                storage.mode(y) <- "double"
                interpx <- rep(0, 2 * n)
                storage.mode(interpx) <- "double"
                interpy <- rep(0, 2 * n)
                storage.mode(interpy) <- "double"
                datatyp <- matrix(0, n, 3)
                storage.mode(datatyp) <- "double"
                datatyp2 <- matrix(0, n, 2)
                storage.mode(datatyp2) <- "double"
                pxpy <- matrix(0, n, 3)
                storage.mode(pxpy) <- "double"
                whisk <- 2
                abagplot.uit <- .Fortran("abagplot", as.integer(n), 
                  as.integer(alpha1), x, y, as.integer(whisk), 
                  tukm = double(2), interpx = interpx, interpy = interpy, 
                  num = as.integer(0), datatyp = datatyp, indoutl = integer(n), 
                  datatyp2 = datatyp2, pxpy = pxpy, boxpl = as.integer(0), 
                  nointer = as.integer(0), PACKAGE = "UBbipl")
                tukmedian <- abagplot.uit$tukm
                x.vec <- abagplot.uit$interpx
                y.vec <- abagplot.uit$interpy
                if (all(x.vec == 0) & all(y.vec == 0)) 
                  stop(message = " x and y both null vectors")
                nie.nul <- !((x.vec == 0) & (y.vec == 0))
                if (Tukey.median) {
                  points(x = tukmedian[1], y = tukmedian[2], 
                    pch = 10, cex = 0.9, col = as.vector(colr))
                }
                polygon(x.vec[nie.nul], y.vec[nie.nul], density = 0, 
                  col = colr, lwd = line.width, lty = as.vector(lty))
            }
        }
        return(legend.labs.bags)
    }
    if (!all(factor.type %in% c("nom", "ord"))) 
        stop("factor.type is a vector of 'nom' and 'ord' \n")
    n <- nrow(Xcat)
    p <- ncol(Xcat)
    pos <- pos[1]
    Gk <- vector("list", p)
    Lk.mat <- vector("list", p)
    Lk <- rep(NA, p)
    zlist <- vector("list", p)
    Hmat <- matrix(nrow = n, ncol = p)
    CatLabels <- NULL
    for (k in 1:p) {
        Gk[[k]] <- indmat(Xcat[, k])
        CatLabels <- append(CatLabels, dimnames(Gk[[k]])[[2]])
        Lk[k] <- ncol(Gk[[k]])
        Lk.mat[[k]] <- t(Gk[[k]]) %*% Gk[[k]]
    }
    G <- matrix(unlist(Gk), nrow = n)
    L.mat <- diag(diag(t(G) %*% G))
    L.mat.min.half <- diag(diag(L.mat)^-0.5)
    N.mat <- diag(n) - matrix(1, nrow = n, ncol = n)/n
    svd.mca <- svd(L.mat.min.half %*% t(G) %*% N.mat %*% G %*% 
        L.mat.min.half/p)
    z.ini <- svd.mca$v[, 1]
    names(z.ini) <- CatLabels
    count <- 0
    for (k in 1:p) {
        zlist[[k]] <- as.vector(z.ini[(count + 1):(count + ncol(Gk[[k]]))])
        names(zlist[[k]]) <- CatLabels[(count + 1):(count + ncol(Gk[[k]]))]
        count <- count + ncol(Gk[[k]])
        plot(1:length(zlist[[k]]), zlist[[k]], type = "l")
        if (identical(factor.type[k], "ord")) {
            zlist[[k]] <- ties(d = 1:length(zlist[[k]]), delta = zlist[[k]])
            lines(1:length(zlist[[k]]), zlist[[k]], col = "red")
        }
        if (identical(factor.type[k], "nom")) {
            lines(1:length(zlist[[k]]), zlist[[k]], col = "green")
        }
        dev.new()
        par(pty = "s", mar = parplotmar)
    }
    for (k in 1:p) {
        labs <- names(zlist[[k]])
        if (identical(factor.type[k], "nom")) 
            zlist[[k]] <- as.vector(zlist[[k]]/sqrt(t(zlist[[k]]) %*% 
                t(Gk[[k]]) %*% N.mat %*% Gk[[k]] %*% zlist[[k]]))
        if (identical(factor.type[k], "ord")) {
            u.temp <- sum(diag(Lk.mat[[k]]))
            v.temp <- sum(zlist[[k]] * diag(Lk.mat[[k]]))
            w.temp <- t(zlist[[k]]) %*% Lk.mat[[k]] %*% zlist[[k]]
            b.temp <- sqrt(1/(w.temp - v.temp^2/u.temp))
            zlist[[k]] <- b.temp * zlist[[k]] - b.temp * v.temp/u.temp
        }
        names(zlist[[k]]) <- labs
        Hmat[, k] <- Gk[[k]] %*% zlist[[k]]
    }
    Hmat <- as.matrix(Hmat)
    dimnames(Hmat) <- dimnames(Xcat)
    again <- TRUE
    RSS.old <- NULL
    Hmat.cent <- N.mat %*% Hmat
    iter <- 0
    while (again) {
        iter <- iter + 1
        svd.out <- svd(Hmat.cent)
        Yr <- scale(svd.out$u[, 1:r] %*% diag(svd.out$d[1:r]) %*% 
            t(svd.out$v[, 1:r]), scale = FALSE, center = TRUE)
        RSS <- sum((Hmat.cent - Yr)^2)
        if (!is.null(RSS.old)) 
            again <- ((RSS.old - RSS) > epsilon)
        RSS.old <- RSS
        if (again) {
            for (k in 1:p) {
                labs <- names(zlist[[k]])
                if (identical(factor.type[k], "nom")) {
                  scale.factor <- as.numeric(sqrt(t(Yr[, k]) %*% 
                    Gk[[k]] %*% solve(Lk.mat[[k]]) %*% t(Gk[[k]]) %*% 
                    Yr[, k]))
                  zlist[[k]] <- as.vector((solve(Lk.mat[[k]]) %*% 
                    t(Gk[[k]]) %*% Yr[, k])/scale.factor)
                }
                if (identical(factor.type[k], "ord")) {
                  zlist[[k]] <- as.vector((solve(Lk.mat[[k]]) %*% 
                    t(Gk[[k]]) %*% Yr[, k]))
                  zlist[[k]] <- ties(d = 1:length(zlist[[k]]), 
                    delta = zlist[[k]])
                  u.temp <- sum(diag(Lk.mat[[k]]))
                  v.temp <- sum(zlist[[k]] * diag(Lk.mat[[k]]))
                  w.temp <- t(zlist[[k]]) %*% Lk.mat[[k]] %*% 
                    zlist[[k]]
                  b.temp <- sqrt(1/(w.temp - v.temp^2/u.temp))
                  zlist[[k]] <- b.temp * zlist[[k]] - b.temp * 
                    v.temp/u.temp
                }
                names(zlist[[k]]) <- labs
                Hmat[, k] <- Gk[[k]] %*% zlist[[k]]
            }
        }
        Hmat.cent <- N.mat %*% Hmat
    }
    Lz <- rep(NA, p)
    for (i in 1:p) Lz[i] <- sum(Lk.mat[[i]] %*% zlist[[i]])
    zLz <- rep(NA, p)
    for (i in 1:p) zLz[i] <- sum(t(zlist[[i]]) %*% Lk.mat[[i]] %*% 
        zlist[[i]])
    Hmat.cent.svd <- svd(Hmat.cent)
    par(pty = "s", mar = parplotmar)
    vars.points <- Hmat.cent.svd$v[, 1:2]
    samples.points <- Hmat.cent %*% Hmat.cent.svd$v[, 1:2]
    grad <- rep(NA, p)
    for (i in 1:p) {
        regr.line <- lm(c(0, vars.points[i, 2]) ~ c(0, vars.points[i, 
            1]))
        grad[i] <- coefficients(regr.line)[2]
    }
    means <- apply(Hmat.cent, 2, mean)
    sd <- sqrt(apply(Hmat.cent, 2, var))
    axes.rows <- 1/(diag(vars.points %*% t(vars.points))) * vars.points
    z.axes <- lapply(1:p, function(j, X.1, axes.rows, orthog.transx, 
        orthog.transy) {
        phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
            axes.rows %*% c(orthog.transx[j], orthog.transy[j])
        axis.vals <- (X.1[, j])
        axis.points <- matrix(0, nrow = length(axis.vals), ncol = 4)
        axis.points[, 1] <- orthog.transx[j] + (axis.vals - phi.vec[j]) * 
            axes.rows[j, 1]
        axis.points[, 2] <- orthog.transy[j] + (axis.vals - phi.vec[j]) * 
            axes.rows[j, 2]
        axis.points[, 3] <- axis.vals
        axis.points[, 4] <- NA
        axis.points <- axis.points[order(axis.points[, 1]), ]
        return(unique(zapsmall(axis.points)))
    }, X.1 = Hmat.cent, axes.rows = axes.rows, orthog.transx = orthog.transx, 
        orthog.transy = orthog.transy)
    all.axes.points <- rbind(z.axes[[1]][, 1:2], z.axes[[2]][, 
        1:2])
    for (i in 3:p) all.axes.points <- rbind(all.axes.points, 
        z.axes[[i]][, 1:2])
    all.points <- samples.points
    plot(all.points[, 1] * exp.factor, all.points[, 2] * exp.factor, 
        xlim = range(all.points[, 1] * exp.factor), ylim = range(all.points[, 
            2] * exp.factor), xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", type = "n", xaxs = "i", yaxs = "i", asp = 1, 
        bty = boxtype)
    if (is.null(ax)) 
        axes <- NULL
    else axes <- (1:p)[ax]
    for (i in axes) {
        z.axes[[i]] <- as.data.frame(z.axes[[i]])
        if (factor.type[i] == "nom") 
            z.axes[[i]][, 4] <- names(zlist[[i]])[match(round(z.axes[[i]][, 
                3], digits = 6), round(zlist[[i]], digits = 6))]
        if (factor.type[i] == "ord") {
            labs.temp <- tied.names(names(zlist[[i]]), zlist[[i]])
            zlist.temp <- unique(zlist[[i]])
            names(zlist.temp) <- labs.temp
            z.axes[[i]][, 4] <- names(zlist.temp)[match(round(z.axes[[i]][, 
                3], digits = 6), round(zlist.temp, digits = 6))]
            labs.temp <- NULL
            zlist.temp <- NULL
        }
    }
    for (i in axes) {
        Draw.line2(x.vals = z.axes[[i]][, 1], y.vals = z.axes[[i]][, 
            2], marker.vals = z.axes[[i]][, 3], line.name = dimnames(Xcat)[[2]][i], 
            axis.col = "white", pos = pos, ax.name.size = ax.name.size, 
            offset = offset)
        n.cat <- nrow(z.axes[[i]])
        if (factor.type[i] == "nom") 
            draw.halfway.lines(x = z.axes[[i]][, 1], y = z.axes[[i]][, 
                2], colour = nom.col[1:n.cat], width = rep(w.factor, 
                n.cat))
        if (factor.type[i] == "ord") {
            temp <- z.axes[[i]][, 1:2]
            draw.halfway.lines(x = temp[, 1], y = temp[, 2], 
                colour = ord.col[1:n], width = 1:(n.cat + 1), 
                w.factor = w.factor, reverse = reverse)
        }
        points(z.axes[[i]][, 1:2], pch = calibration.pch, col = calibration.col, 
            cex = calibration.size)
        for (j in 1:n.cat) {
            text(x = z.axes[[i]][j, 1], y = z.axes[[i]][j, 2], 
                label = z.axes[[i]][j, 4], cex = calibration.label.size, 
                col = calibration.label.col, pos = calibration.label.pos[i], 
                offset = calibration.label.offset[i], xpd = NA)
        }
    }
    if (!is.null(plot.samples)) {
        points(samples.points[plot.samples, ], pch = samples.pch, 
            col = samples.col, cex = samples.size)
        if (label.samples) {
            sample.labels <- ID[plot.samples]
            text(samples.points[plot.samples, ], labels = sample.labels, 
                pos = sample.label.pos, cex = sample.label.size, 
                offset = sample.label.offset)
        }
    }
    if (!is.null(drawbagplots)) {
        usr <- par("usr")
        .bags.plot(Z = samples.points, specify.bags = specify.bags, 
            line.type.bags = line.type.bags, class.vec = class.vec, 
            class.pch = class.pch, class.cols = class.cols, usr = usr, 
            alpha1 = alpha, max.num = max.num, line.width = line.width, 
            Tukey.median = FALSE, c.hull.n = c.hull.n)
    }
    if (!is.null(predict.sample)) {
        for (j in predict.sample) {
            for (i in 1:p) DrawOrthogline(x1 = z.axes[[i]][1, 
                1], x2 = z.axes[[i]][2, 1], y1 = z.axes[[i]][1, 
                2], y2 = z.axes[[i]][2, 2], px = samples.points[j, 
                1], py = samples.points[j, 2], ort.lty = ort.lty)
        }
    }
    if (select.origin) {
        cat("Move cursor where you want origin and press left mouse button \n")
        flush.console()
        origin.pos <- locator(1)
        z.axes <- lapply(1:p, function(j, X.1, axes.rows, orthog.transx, 
            orthog.transy) {
            phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                axes.rows %*% c(orthog.transx[j], orthog.transy[j])
            axis.vals <- (X.1[, j])
            axis.points <- matrix(0, nrow = length(axis.vals), 
                ncol = 4)
            axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 1]
            axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 2]
            axis.points[, 3] <- axis.vals
            axis.points[, 4] <- NA
            axis.points <- axis.points[order(axis.points[, 1]), 
                ]
            return(unique(zapsmall(axis.points)))
        }, X.1 = Hmat.cent, axes.rows = axes.rows, orthog.transx = rep(origin.pos$x, 
            p), orthog.transy = rep(origin.pos$y, p))
        all.axes.points <- rbind(z.axes[[1]][, 1:2], z.axes[[2]][, 
            1:2])
        for (i in 3:p) all.axes.points <- rbind(all.axes.points, 
            z.axes[[i]][, 1:2])
        all.points <- rbind(all.axes.points, samples.points)
        plot(all.points[, 1] * exp.factor, all.points[, 2] * 
            exp.factor, xlim = range(all.points[, 1] * exp.factor), 
            ylim = range(all.points[, 2] * exp.factor), xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", type = "n", xaxs = "i", 
            yaxs = "i", asp = 1)
        if (is.null(ax)) 
            axes <- NULL
        else axes <- (1:p)[ax]
        for (i in axes) {
            z.axes[[i]] <- as.data.frame(z.axes[[i]])
            if (factor.type[i] == "nom") 
                z.axes[[i]][, 4] <- names(zlist[[i]])[match(round(z.axes[[i]][, 
                  3], digits = 6), round(zlist[[i]], digits = 6))]
            if (factor.type[i] == "ord") {
                labs.temp <- tied.names(names(zlist[[i]]), zlist[[i]])
                zlist.temp <- unique(zlist[[i]])
                names(zlist.temp) <- labs.temp
                z.axes[[i]][, 4] <- names(zlist.temp)[match(round(z.axes[[i]][, 
                  3], digits = 6), round(zlist.temp, digits = 6))]
                labs.temp <- NULL
                Zlist.temp <- NULL
            }
        }
        for (i in axes) {
            Draw.line2(x.vals = z.axes[[i]][, 1], y.vals = z.axes[[i]][, 
                2], marker.vals = z.axes[[i]][, 3], line.name = dimnames(Xcat)[[2]][i], 
                axis.col = "white", ax.name.size = ax.name.size, 
                pos = pos, offset = offset)
            n.cat <- nrow(z.axes[[i]])
            if (factor.type[i] == "nom") 
                draw.halfway.lines(x = z.axes[[i]][, 1], y = z.axes[[i]][, 
                  2], colour = nom.col[1:n.cat], width = rep(w.factor, 
                  n.cat))
            if (factor.type[i] == "ord") {
                temp <- z.axes[[i]][, 1:2]
                draw.halfway.lines(x = temp[, 1], y = temp[, 
                  2], colour = ord.col[1:n], width = 1:(n.cat + 
                  1), w.factor = w.factor, reverse = reverse)
            }
            points(z.axes[[i]][, 1:2], pch = calibration.pch, 
                col = calibration.col, cex = calibration.size)
            for (j in 1:n.cat) {
                text(x = z.axes[[i]][j, 1], y = z.axes[[i]][j, 
                  2], label = z.axes[[i]][j, 4], cex = calibration.label.size, 
                  col = calibration.label.col, pos = calibration.label.pos[i], 
                  offset = calibration.label.offset[i], xpd = NA)
            }
        }
        if (!is.null(plot.samples)) {
            points(samples.points[plot.samples, ], pch = samples.pch, 
                col = samples.col, cex = samples.size)
            if (label.samples) {
                sample.labels <- ID[plot.samples]
                text(samples.points[plot.samples, ], labels = sample.labels, 
                  pos = sample.label.pos, cex = sample.label.size, 
                  offset = sample.label.offset)
            }
        }
        if (!is.null(predict.sample)) {
            for (j in predict.sample) {
                for (i in 1:p) DrawOrthogline(x1 = z.axes[[i]][1, 
                  1], x2 = z.axes[[i]][2, 1], y1 = z.axes[[i]][1, 
                  2], y2 = z.axes[[i]][2, 2], px = samples.points[j, 
                  1], py = samples.points[j, 2], ort.lty = ort.lty)
            }
        }
    }
    points(0, 0, pch = 3, cex = 2)
    dev.new()
    PCAbipl.cat(Hmat.cent, exp.factor = exp.factor)
    dev.new()
    par(mfrow = z.score.graph)
    for (k in 1:p) {
        par(pty = "m", mar = parplotmar)
        plot(1:length(zlist[[k]]), zlist[[k]], type = "b", xaxt = "n", 
            ylab = "Final z quantification", xlab = "", main = dimnames(Xcat)[[2]][k])
        axis(side = 1, at = 1:length(names(zlist[[k]])), labels = names(zlist[[k]]))
    }
    par(mfrow = c(1, 1))
    list(iter = iter, Hmat = Hmat, Hmat.cent = Hmat.cent, Yr = Yr, 
        zlist = zlist, z.axes = z.axes, samples.points = samples.points)
}
