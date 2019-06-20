drawbipl.bagalpha.density.groups <-
function (Z, z.axes, z.axes.names = NULL, alpha = 0.5, basket.beta = 0.5, 
    basket.n = 36, ax = NULL, ax.col, ax.name.col, ax.name.size, 
    bandwidth, pch.samples.size = 1, c.hull.n = 10, class.vec, 
    colours.density, cuts.density, dens.ax.cex, dens.ax.tcl, 
    dens.mgp, draw.densitycontours, exp.factor = 1.2, label = TRUE, 
    label.size, large.scale = FALSE, layout.heights, line.width = 1, 
    marker.size, markers = TRUE, max.num, means.plot = FALSE, 
    offset = rep(0.5, 4), offset.m, p, parlegendmar.dens, pch.means = 15, 
    pch.means.size = 1, pos, pos.m, side.label = side.label, 
    smooth.n, specify.bags, specify.classes, specify.density.classes, 
    specify.beta.baskets, strepie = c(1, 1), Title = NULL, Tukey.median = TRUE, 
    Z.means.mat) 
{
    .axes.plot <- function(z.axes, z.axes.names, p, ax, ax.col, 
        ax.name.col, ax.name.size, markers, usr, label, label.size, 
        marker.size, offset, offset.m, pos, pos.m, side.label, 
        strepie) {
        if (is.null(ax)) 
            axes <- NULL
        else axes <- (1:p)[ax]
        for (i in axes) {
            marker.mat <- z.axes[[i]][z.axes[[i]][, 4] == 1, 
                1:3]
            x.vals <- marker.mat[, 1]
            y.vals <- marker.mat[, 2]
            if (is.null(z.axes.names)) {
                axis.name <- paste("v", i, sep = "")
            }
            else {
                axis.name <- z.axes.names[i]
            }
            std.markers <- marker.mat[, 3]
            x.invals <- x.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            if (length(x.invals) < 2) 
                warning(paste("Less than 2 markers on axis ", 
                  i, ". Increase n.int."))
            y.invals <- y.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            tick.labels <- zapsmall(std.markers[x.vals < usr[2] & 
                x.vals > usr[1] & y.vals < usr[4] & y.vals > 
                usr[3]])
            uit <- Draw.line2(x.vals = x.invals, y.vals = y.invals, 
                marker.vals = tick.labels, line.name = axis.name, 
                ax.name.size = ax.name.size, axis.col = ax.col$ax.col[i], 
                ax.name.col = ax.name.col[i], offset = offset, 
                pos = pos)
            gradient <- uit$gradient
            for (j in 1:length(x.invals)) Draw.onecmline(x = x.invals[j], 
                y = y.invals[j], grad = -1/gradient, expand = strepie[1], 
                both.sides = TRUE, col = ax.col$tickmarker.col[i])
            if (markers == TRUE) {
                x.labvals <- x.invals
                y.labvals <- y.invals
            }
            else {
                x.labvals <- x.invals[c(1, length(x.invals))]
                y.labvals <- y.invals[c(1, length(y.invals))]
                tick.labels <- tick.labels[c(1, length(tick.labels))]
            }
            for (j in 1:length(x.labvals)) Plot.marker.new(x = x.labvals[j], 
                y = y.labvals[j], grad = -1/gradient, mark = tick.labels[j], 
                expand = strepie[2], marker.size = marker.size, 
                col = ax.col$marker.col[i], offset.m = offset.m[i], 
                pos.m = NULL, side.label = side.label[i])
        }
    }
    .samples.plot <- function(Z, pch.samples.size, pch.means, 
        Z.means.mat, specify.classes, usr, label, label.size, 
        class.vec) {
        if (ncol(Z) == 2) 
            Z <- data.frame(Z, pch.means = pch.means, stringsAsFactors = FALSE)
        if (ncol(Z) == 3) 
            Z <- data.frame(Z, colr = 1, stringsAsFactors = FALSE)
        if (ncol(Z) == 4) 
            Z <- data.frame(Z, line.type = 1, stringsAsFactors = FALSE)
        if (ncol(Z) == 5) 
            Z <- cbind(Z, pch.samples.size = pch.samples.size, 
                stringsAsFactors = FALSE)
        Z.plot <- Z
        x.vals <- Z.plot[, 1]
        y.vals <- Z.plot[, 2]
        invals <- x.vals < usr[2] & x.vals > usr[1] & y.vals < 
            usr[4] & y.vals > usr[3]
        Z.plot <- Z.plot[invals, ]
        if (is.null(Z.means.mat)) 
            J <- 0
        else J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[specify.classes]
        legend.labs <- classes
        for (j in classes) {
            Z.class <- data.frame(Z[class.vec == j, , drop = FALSE], 
                stringsAsFactors = FALSE)
            if (label == TRUE) 
                text(Z.class[, 1], Z.class[, 2] - 0.015 * (usr[4] - 
                  usr[3]), labels = dimnames(Z.class)[[1]], cex = label.size)
            chr <- Z.class[1, 3]
            colr <- Z.class[1, 4]
            lty <- Z.class[1, 5]
            if (!is.null(specify.classes)) 
                for (i in 1:nrow(Z.class)) points(x = Z.class[i, 
                  1], y = Z.class[i, 2], pch = Z.class[i, 3], 
                  col = Z.class[i, 4], lty = Z.class[i, 5], cex = Z.class[i, 
                    6])
        }
    }
    .bags.plot <- function(Z, Z.means.mat, specify.bags, class.vec, 
        usr, alpha1, max.num, line.width, Tukey.median, c.hull.n) {
        J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[specify.bags]
        legend.labs.bags <- classes
        for (j in classes) {
            Z.class <- as.data.frame(Z[class.vec == j, , drop = FALSE], 
                stringsAsFactors = FALSE)
            chr <- Z.class[1, 3]
            colr <- Z.class[1, 4]
            lty <- Z.class[1, 5]
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
    .beta.basket.plot <- function(Z, Z.means.mat, specify.beta.baskets, 
        class.vec, line.width, basket.beta, colr, basket.n) {
        J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[specify.beta.baskets]
        legend.labs.baskets <- classes
        for (j in classes) {
            Z.class <- Z[class.vec == j, ]
            chr <- Z.class[1, 3]
            colr <- as.character(Z.class[1, 4])
            lty <- Z.class[1, 5]
            flush.console()
            cat(paste("basket for class ", j, " with ", nrow(Z.class), 
                " samples in class ", j, sep = ""), "\n")
            x <- Z.class[, 1]
            y <- Z.class[, 2]
            ScatterplotBaskets(x = x, y = y, n = basket.n, pp = basket.beta, 
                lty = lty, col = colr, lwd = line.width)
        }
        return(legend.labs.baskets)
    }
    if (alpha < 0 | alpha > 0.99) 
        stop(message = "alpha not to be negative or larger than 0.99")
    alpha.entered <- alpha
    alpha <- round(alpha, digits = 2)
    if (abs(alpha.entered - alpha) > 0) 
        cat("alpha has been rounded to ", alpha, "\n")
    p <- length(z.axes)
    alpha1 <- 100 * alpha
    if (means.plot == FALSE & is.null(specify.bags) & is.null(specify.classes)) 
        stop("You have specified nothing to be plotted")
    tel <- 1
    if (large.scale) {
        specify.bags <- NULL
        specify.classes <- NULL
        legend.lab.bags <- NULL
        specify.beta.baskets <- NULL
        density.classes <- unique(dimnames(Z.means.mat)[[1]])[specify.density.classes]
        Z.plot <- data.frame(Z.means.mat, pch.means.size = pch.means.size)
        for (i in density.classes) {
            dev.new()
            Z.class <- Z[class.vec == i, ]
            layout(matrix(1:2, ncol = 1), heights = layout.heights)
            par(pty = "s", cex = 1)
            plot(Z.plot[, 1], Z.plot[, 2], xlim = range(Z.plot[, 
                1]) * exp.factor, ylim = range(Z.plot[, 2]) * 
                exp.factor, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", type = "n", xaxs = "i", yaxs = "i", 
                asp = 1)
            usr <- par("usr")
            par(new = TRUE, pty = "s")
            if (!is.null(bandwidth)) 
                ff1 <- kde2d(Z.class[, 1], Z.class[, 2], n = smooth.n, 
                  lims = usr, h = bandwidth)
            else ff1 <- kde2d(Z.class[, 1], Z.class[, 2], n = smooth.n, 
                lims = usr)
            levels.rect <- pretty(range(ff1$z), n = cuts.density)
            col.use <- colorRampPalette(colours.density[tel, 
                ])
            tel <- tel + 1
            col.use <- col.use(length(levels.rect) - 1)
            image(ff1, xlim = usr[1:2], ylim = usr[3:4], breaks = levels.rect, 
                xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
                xlab = "", ylab = "", asp = 1, col = col.use)
            par(new = TRUE, pty = "s")
            plot(Z.plot[, 1:2], xlim = usr[1:2], ylim = usr[3:4], 
                xaxs = "i", yaxs = "i", asp = 1, xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", type = "n")
            if (means.plot == TRUE) {
                for (i in 1:nrow(Z.plot)) points(x = Z.plot[i, 
                  1], y = Z.plot[i, 2], pch = Z.plot[i, 3], col = Z.plot[i, 
                  4], cex = Z.plot[i, 6])
                legend.labs <- unique(dimnames(Z.means.mat)[[1]])
                .axes.plot(z.axes = z.axes, z.axes.names = z.axes.names, 
                  p = p, ax = ax, ax.name.col = ax.name.col, 
                  ax.name.size = ax.name.size, ax.col = ax.col, 
                  markers = markers, usr = usr, label = label, 
                  label.size = label.size, marker.size = marker.size, 
                  offset = offset, offset.m = offset.m, pos = pos, 
                  pos.m = pos.m, side.label = side.label, strepie = strepie)
                warning("When means are plotted on a large scale graph no bags or samples are plotted\n")
            }
            title(main = (ifelse(is.null(Title), paste(alpha1, 
                " % BAG PLOT", sep = ""), Title)))
            par(pty = "m", mar = parlegendmar.dens)
            plot(range(levels.rect), y = 1:2, ylim = c(10, 100), 
                xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", type = "n", frame.plot = FALSE)
            rect(xleft = levels.rect[-length(levels.rect)], ybottom = 10, 
                xright = levels.rect[-1], ytop = 50, col = col.use, 
                border = FALSE)
            axis(side = 1, at = pretty(levels.rect, n = 8), labels = pretty(levels.rect, 
                n = 8), line = 0, cex.axis = dens.ax.cex, mgp = dens.mgp, 
                tcl = dens.ax.tcl, las = 0)
        }
    }
    else {
        density.classes <- unique(dimnames(Z.means.mat)[[1]])[specify.density.classes]
        for (i in density.classes) {
            dev.new()
            Z.class <- Z[class.vec == i, ]
            layout(matrix(1:2, ncol = 1), heights = layout.heights)
            par(pty = "s", cex = 1)
            plot(x = Z[, 1], y = Z[, 2], xlim = range(Z[, 1]) * 
                exp.factor, ylim = range(Z[, 2]) * exp.factor, 
                xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
                type = "n", xaxs = "i", yaxs = "i", asp = 1)
            usr <- par("usr")
            par(new = TRUE, pty = "s")
            if (!is.null(bandwidth)) 
                ff1 <- kde2d(Z.class[, 1], Z.class[, 2], n = smooth.n, 
                  lims = usr, h = bandwidth)
            else ff1 <- kde2d(Z.class[, 1], Z.class[, 2], n = smooth.n, 
                lims = usr)
            levels.rect <- pretty(range(ff1$z), n = cuts.density)
            col.use <- colorRampPalette(colours.density[tel, 
                ])
            tel <- tel + 1
            col.use <- col.use(length(levels.rect) - 1)
            image(ff1, xlim = usr[1:2], ylim = usr[3:4], breaks = levels.rect, 
                xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
                xlab = "", ylab = "", asp = 1, col = col.use)
            par(new = TRUE, pty = "s")
            plot(Z[, 1:2], xlim = usr[1:2], ylim = usr[3:4], 
                xaxs = "i", yaxs = "i", asp = 1, xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", type = "n")
            legend.lab.bags <- NULL
            if (means.plot == TRUE & large.scale == TRUE) {
                for (i in 1:nrow(Z.plot)) points(x = Z.plot[i, 
                  1], y = Z.plot[i, 2], pch = Z.plot[i, 3], col = Z.plot[i, 
                  4], cex = Z.plot[i, 6])
                legend.labs <- unique(dimnames(Z.means.mat)[[1]])
                .axes.plot(z.axes = z.axes, z.axes.names = z.axes.names, 
                  p = p, ax = ax, ax.name.col = ax.name.col, 
                  ax.name.size = ax.name.size, ax.col = ax.col, 
                  markers = markers, usr = usr, label = label, 
                  label.size = label.size, marker.size = marker.size, 
                  offset = offset, offset.m = offset.m, pos = pos, 
                  pos.m = pos.m, side.label = side.label, strepie = strepie)
                warning("When means are plotted on a large scale graph no bags or samples are plotted\n")
            }
            if (means.plot == TRUE & large.scale == FALSE) {
                Z.plot <- cbind(Z.means.mat, pch.means.size)
                for (i in 1:nrow(Z.plot)) points(x = Z.plot[i, 
                  1], y = Z.plot[i, 2], pch = Z.plot[i, 3], col = Z.plot[i, 
                  4], cex = Z.plot[i, 6])
                legend.labs <- unique(dimnames(Z.means.mat)[[1]])
                .axes.plot(z.axes = z.axes, z.axes.names = z.axes.names, 
                  p = p, ax = ax, ax.name.col = ax.name.col, 
                  ax.name.size = ax.name.size, ax.col = ax.col, 
                  markers = markers, usr = usr, label = label, 
                  label.size = label.size, marker.size = marker.size, 
                  offset = offset, offset.m = offset.m, pos = pos, 
                  pos.m = pos.m, side.label = side.label, strepie = strepie)
            }
            if (means.plot == FALSE & large.scale == FALSE) {
                legend.labs <- unique(dimnames(Z.means.mat)[[1]])
                .axes.plot(z.axes = z.axes, z.axes.names = z.axes.names, 
                  p = p, ax = ax, ax.name.col = ax.name.col, 
                  ax.name.size = ax.name.size, ax.col = ax.col, 
                  markers = markers, usr = usr, label = label, 
                  label.size = label.size, marker.size = marker.size, 
                  offset = offset, offset.m = offset.m, pos = pos, 
                  pos.m = pos.m, side.label = side.label, strepie = strepie)
            }
            if (!is.null(specify.classes)) 
                .samples.plot(Z = Z, pch.samples.size = pch.samples.size, 
                  pch.means = pch.means, Z.means.mat = Z.means.mat, 
                  specify.classes = specify.classes, usr = usr, 
                  label = label, label.size = label.size, class.vec = class.vec)
            if (!is.null(specify.bags)) 
                legend.lab.bags <- .bags.plot(Z = Z, Z.means.mat = Z.means.mat, 
                  specify.bags = specify.bags, usr = usr, class.vec = class.vec, 
                  alpha1 = alpha1, max.num = max.num, line.width = line.width, 
                  Tukey.median = Tukey.median, c.hull.n = c.hull.n)
            if (!is.null(specify.beta.baskets)) 
                legend.lab.beta.baskets <- .beta.basket.plot(Z = Z, 
                  Z.means.mat = Z.means.mat, specify.beta.baskets = specify.beta.baskets, 
                  class.vec = class.vec, basket.beta = basket.beta, 
                  basket.n = basket.n, line.width = line.width, 
                  col = colr)
            title(main = (ifelse(is.null(Title), paste(alpha1, 
                " % BAG PLOT", sep = ""), Title)))
        }
    }
    list(legend.lab.bags = legend.lab.bags, col.use = col.use, 
        levels.rect.rect = levels.rect)
}
