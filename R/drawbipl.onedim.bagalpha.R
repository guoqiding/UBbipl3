drawbipl.onedim.bagalpha <-
function (Z, z.axes, z.axes.names = NULL, ax = NULL, ax.col, 
    ax.name.col, ax.name.size, alpha = 0.5, c.hull.n = 10, class.vec, 
    constant = 0, conf.alpha = NULL, density.plot = FALSE, exp.factor = 1.2, 
    label = TRUE, label.size, large.scale = FALSE, line.width = 1, 
    markers = TRUE, marker.size, max.num, means.plot = FALSE, 
    offset = rep(0.5, 4), offset.m, ort.lty = 1, p, pch.means = 15, 
    pch.means.size = 1, pch.samples.size = 1, pos, pos.m, predictions.mean = NULL, 
    predictions.sample = NULL, side.label = side.label, specify.bags, 
    specify.classes, strepie = c(1, 1), Title = NULL, Tukey.median = TRUE, 
    Z.means.mat) 
{
    .axes.plot <- function(Z, Z.means.mat, z.axes, z.axes.names, 
        p, ax, ax.col, ax.name.col, ax.name.size, markers, usr, 
        label, label.size, marker.size, offset, offset.m, pos, 
        pos.m, side.label, strepie, predictions.sample, predictions.mean, 
        constant, ort.lty) {
        if (is.null(ax)) 
            axes <- NULL
        else axes <- (1:p)[ax]
        if (!is.null(predictions.sample)) {
            predictions <- data.frame(matrix(NA, nrow = length(axes), 
                ncol = length(predictions.sample)))
            dimnames(predictions) <- list(1:nrow(predictions), 
                paste("s", predictions.sample, sep = ""))
        }
        if (!is.null(predictions.mean)) {
            predictions <- data.frame(matrix(NA, nrow = length(axes), 
                ncol = length(predictions.mean)))
            dimnames(predictions) <- list(1:nrow(predictions), 
                paste("m", predictions.mean, sep = ""))
        }
        r.names = NULL
        for (i in axes) {
            z.axes[[i]][, 2] <- constant * i
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
            r.names[i] <- axis.name
            std.markers <- zapsmall(marker.mat[, 3])
            x.invals <- x.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            if (length(x.invals) < 2) {
                warning(paste("Less than 2 markers on axis ", 
                  i, ". Increase n.int or check value of .constant."))
                marker.mat <- z.axes[[i]][, 1:3]
                x.vals <- marker.mat[, 1]
                y.vals <- marker.mat[, 2]
                std.markers <- zapsmall(marker.mat[, 3])
                x.invals <- x.vals[x.vals < usr[2] & x.vals > 
                  usr[1] & y.vals < usr[4] & y.vals > usr[3]]
                y.invals <- y.vals[x.vals < usr[2] & x.vals > 
                  usr[1] & y.vals < usr[4] & y.vals > usr[3]]
                tick.labels <- std.markers[x.vals < usr[2] & 
                  x.vals > usr[1] & y.vals < usr[4] & y.vals > 
                  usr[3]]
                Draw.line2(x.vals = x.invals, y.vals = y.invals, 
                  marker.vals = tick.labels, line.name = axis.name, 
                  ax.name.size = ax.name.size, axis.col = ax.col$ax.col[i], 
                  ax.name.col = ax.name.col[i], offset = offset, 
                  pos = pos)
                next
            }
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
            if (!is.null(predictions.sample)) {
                for (ss in predictions.sample) {
                  predictions[i, paste("s", ss, sep = "")] <- round(DrawOrthogline(x1 = x.invals[1], 
                    y1 = y.invals[1], x2 = x.invals[length(x.invals)], 
                    y2 = y.invals[length(y.invals)], val1 = tick.labels[1], 
                    val2 = tick.labels[length(x.invals)], px = Z[ss, 
                      1], py = Z[ss, 2], ort.lty = ort.lty), 
                    digits = 4)
                }
            }
            if (!is.null(predictions.mean)) {
                for (ss in predictions.mean) {
                  predictions[i, paste("m", ss, sep = "")] <- round(DrawOrthogline(x1 = x.invals[1], 
                    y1 = y.invals[1], x2 = x.invals[length(x.invals)], 
                    y2 = y.invals[length(y.invals)], val1 = tick.labels[1], 
                    val2 = tick.labels[length(x.invals)], px = Z.means.mat[ss, 
                      1], py = Z.means.mat[ss, 2], ort.lty = ort.lty), 
                    digits = 4)
                }
            }
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
            for (j in 1:length(x.labvals)) Plot.marker(x = x.labvals[j], 
                y = y.labvals[j], grad = -1/gradient, mark = tick.labels[j], 
                expand = strepie[2], marker.size = marker.size, 
                col = ax.col$marker.col[i], offset.m = offset.m[i], 
                pos.m = pos.m[i], side.label = side.label[i])
        }
        if (is.null(predictions.sample) & is.null(predictions.mean)) 
            predictions <- NULL
        if (!is.null(predictions)) {
            predictions <- na.omit(predictions)
            predictions <- as.matrix(predictions)
            dimnames(predictions)[[1]] <- r.names[!is.na(r.names)]
        }
        return(predictions)
    }
    .samples.plot <- function(Z, pch.samples.size, pch.means, 
        Z.means.mat, specify.classes, usr, label, label.size, 
        class.vec) {
        if (ncol(Z) == 3) 
            Z <- data.frame(Z, 1, stringsAsFactors = FALSE)
        if (ncol(Z) == 4) 
            Z <- data.frame(Z, 1, stringsAsFactors = FALSE)
        if (ncol(Z) == 5) 
            Z <- data.frame(Z, pch.samples.size, stringsAsFactors = FALSE)
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
                  col = Z.class[i, 4], cex = Z.class[i, 6])
        }
    }
    .bags.plot <- function(Z, Z.means.mat, specify.bags, class.vec, 
        usr, alpha1, max.num, line.width, Tukey.median, c.hull.n) {
        J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[specify.bags]
        legend.labs.bags <- classes
        for (j in classes) {
            Z.class <- data.frame(Z[class.vec == j, ], stringsAsFactors = FALSE)
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
    if (large.scale) {
        if (nrow(Z.means.mat) == 1) 
            stop("This option should only be chosen if there are at least two different classes")
        specify.bags <- NULL
        specify.classes <- specify.classes
        Z.plot <- cbind(Z.means.mat, pch.means.size)
        plot(Z.plot[, 1] * exp.factor, Z.plot[, 2] * exp.factor, 
            xlim = range(Z.plot[, 1] * exp.factor), ylim = range(Z.plot[, 
                2] * exp.factor), xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", type = "n", xaxs = "i", yaxs = "i", asp = 1)
        abline(h = 0, col = "black")
    }
    else plot(Z[, 1] * exp.factor, Z[, 2] * exp.factor, xlim = range(Z[, 
        1] * exp.factor), ylim = range(Z[, 2] * exp.factor), 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", 
        xaxs = "i", yaxs = "i", asp = 1)
    abline(h = 0, col = "black")
    usr <- par("usr")
    legend.lab.bags <- NULL
    J <- nrow(Z.means.mat)
    classes.cc <- rownames(Z.means.mat)
    hh <- 0.1 * (usr[4] - usr[2])
    sign <- sign(hh)
    if (means.plot == TRUE & large.scale == TRUE) {
        if (!is.null(predictions.sample)) {
            predictions.sample <- NULL
            warning("Predictions for sample points only available when large.scale == FALSE \n")
        }
        for (i in 1:nrow(Z.plot)) {
            points(x = Z.plot[i, 1], y = Z.plot[i, 2], pch = Z.plot[i, 
                3], col = Z.plot[i, 4], cex = Z.plot[i, 6])
            if (!is.null(conf.alpha)) {
                Z.class.cc <- Z[class.vec == classes.cc[i], ]
                temp <- sqrt(1/(nrow(Z) - nrow(Z.means.mat))) * 
                  sqrt(1/nrow(Z.class.cc)) * qnorm(conf.alpha + 
                  (1 - conf.alpha)/2)
                x.cc <- Z.means.mat[i, 1]
                y.cc <- Z.means.mat[i, 2]
                lines(x = c(x.cc - temp, x.cc + temp), y = c(y.cc + 
                  2 * hh * sign, y.cc + 2 * hh * sign), col = Z.class.cc[1, 
                  4], lwd = 3)
                lines(x = c(x.cc - temp, x.cc - temp), y = c(y.cc + 
                  hh * sign, y.cc + 3 * hh * sign), col = Z.class.cc[1, 
                  4], lwd = 3)
                lines(x = c(x.cc + temp, x.cc + temp), y = c(y.cc + 
                  hh * sign, y.cc + 3 * hh * sign), col = Z.class.cc[1, 
                  4], lwd = 3)
            }
        }
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
        predictions <- .axes.plot(Z = Z, Z.means.mat = Z.means.mat, 
            z.axes = z.axes, z.axes.names = z.axes.names, p = p, 
            ax = ax, ax.name.size = ax.name.size, ax.col = ax.col, 
            ax.name.col = ax.name.col, markers = markers, usr = usr, 
            label = label, label.size = label.size, marker.size = marker.size, 
            offset = offset, offset.m = offset.m, ort.lty = ort.lty, 
            pos = pos, pos.m = pos.m, side.label = side.label, 
            strepie = strepie, predictions.sample = NULL, predictions.mean = predictions.mean, 
            constant = constant)
        warning("When means are plotted on a large scale graph no bags or samples are plotted\n")
    }
    if (means.plot == TRUE & large.scale == FALSE) {
        if (!is.null(predictions.mean)) {
            warning("Predictions for class means more suitable when large.scale == TRUE \n")
        }
        Z.plot <- cbind(Z.means.mat, pch.means.size)
        for (i in 1:nrow(Z.plot)) {
            points(x = Z.plot[i, 1], y = Z.plot[i, 2], pch = Z.plot[i, 
                3], col = Z.plot[i, 4], cex = Z.plot[i, 6])
            if (!is.null(conf.alpha)) {
                {
                  {
                    Z.class.cc <- Z[class.vec == classes.cc[i], 
                      ]
                    temp <- sqrt(1/(nrow(Z) - nrow(Z.means.mat))) * 
                      sqrt(1/nrow(Z.class.cc)) * qnorm(conf.alpha + 
                      (1 - conf.alpha)/2)
                    x.cc <- Z.means.mat[i, 1]
                    y.cc <- Z.means.mat[i, 2]
                    lines(x = c(x.cc - temp, x.cc + temp), y = c(y.cc + 
                      2 * hh * sign, y.cc + 2 * hh * sign), col = Z.class.cc[1, 
                      4], lwd = 3)
                    lines(x = c(x.cc - temp, x.cc - temp), y = c(y.cc + 
                      hh * sign, y.cc + 3 * hh * sign), col = Z.class.cc[1, 
                      4], lwd = 3)
                    lines(x = c(x.cc + temp, x.cc + temp), y = c(y.cc + 
                      hh * sign, y.cc + 3 * hh * sign), col = Z.class.cc[1, 
                      4], lwd = 3)
                  }
                }
            }
        }
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
        predictions <- .axes.plot(Z = Z, Z.means.mat = Z.means.mat, 
            z.axes = z.axes, z.axes.names = z.axes.names, p = p, 
            ax = ax, ax.name.size = ax.name.size, ax.col = ax.col, 
            ax.name.col = ax.name.col, markers = markers, usr = usr, 
            label = label, label.size = label.size, marker.size = marker.size, 
            offset = offset, offset.m = offset.m, ort.lty = ort.lty, 
            pos = pos, pos.m = pos.m, side.label = side.label, 
            strepie = strepie, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, constant = constant)
    }
    if (means.plot == FALSE & large.scale == FALSE) {
        if (!is.null(predictions.mean)) {
            warning("Predictions for class means only available when means.plot == TRUE \n")
            prediction.means <- NULL
        }
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
        predictions <- .axes.plot(Z = Z, Z.means.mat = Z.means.mat, 
            z.axes = z.axes, z.axes.names = z.axes.names, p = p, 
            ax = ax, ax.name.size = ax.name.size, ax.col = ax.col, 
            ax.name.col = ax.name.col, markers = markers, usr = usr, 
            label = label, label.size = label.size, marker.size = marker.size, 
            offset = offset, offset.m = offset.m, ort.lty = ort.lty, 
            pos = pos, pos.m = pos.m, side.label = side.label, 
            strepie = strepie, predictions.sample = predictions.sample, 
            predictions.mean = NULL, constant = constant)
    }
    if (!is.null(specify.classes)) 
        .samples.plot(Z = Z, pch.samples.size = pch.samples.size, 
            pch.means = pch.means, Z.means.mat = Z.means.mat, 
            specify.classes = specify.classes, usr = usr, label = label, 
            label.size = label.size, class.vec = class.vec)
    if (!is.null(specify.bags)) 
        legend.lab.bags <- .bags.plot(Z = Z, Z.means.mat = Z.means.mat, 
            specify.bags = specify.bags, usr = usr, class.vec = class.vec, 
            alpha1 = alpha1, max.num = max.num, line.width = line.width, 
            Tukey.median = Tukey.median, c.hull.n = c.hull.n)
    title(main = (ifelse(is.null(Title), paste(alpha1, " % BAG PLOT", 
        sep = ""), Title)))
    if (density.plot == "groups") {
        par(pty = "m", new = T, fig = c(0, 1, 0.43, 1))
        names <- unique(dimnames(Z.means.mat)[[1]])[specify.classes]
        Z.points.list <- vector("list", length(names))
        Z.dens.list <- vector("list", length(names))
        max.vec <- rep(0, length(names))
        for (i in 1:length(names)) {
            Z.points.list[[i]] <- Z[class.vec == names[i], ]
            Z.dens.list[[i]] <- density(Z.points.list[[i]][, 
                1])
            max.vec[i] <- max(Z.dens.list[[i]]$y)
        }
        maximum <- max(max.vec)
        plot(Z.dens.list[[1]], xlim = usr[1:2], ylim = c(0, maximum), 
            xaxt = "n", xaxs = "i", xlab = "", ylab = "", main = "", 
            col = Z.points.list[[1]][1, 4], cex.axis = 0.4, lwd = line.width, 
            bty = "n")
        if (length(names) > 1) {
            for (i in 2:length(names)) lines(Z.dens.list[[i]], 
                xlim = usr[1:2], xaxs = "i", col = Z.points.list[[i]][1, 
                  4], lwd = line.width)
        }
    }
    if (density.plot == "all" | density.plot == "TRUE" | density.plot == 
        "T") {
        par(pty = "m", new = T, fig = c(0, 1, 0.43, 1))
        Z.dens.list <- density(Z[, 1])
        max.vec <- max(Z.dens.list$y)
        plot(Z.dens.list, xlim = usr[1:2], ylim = c(0, max.vec), 
            xaxt = "n", xaxs = "i", xlab = "", ylab = "", main = "", 
            col = "black", cex.axis = 0.4, lwd = line.width)
    }
    list(legend.lab.bags = legend.lab.bags, predictions = predictions)
}
