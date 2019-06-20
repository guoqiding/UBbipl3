drawbipl.nonlin <-
function (Z, z.axes, z.axes.names = NULL, p, ax = NULL, ax.col, 
    ax.name.size, alpha = 0.5, pch.samples = 1, pch.samples.size = 1, 
    c.hull.n = 10, class.vec, exp.factor = 1.2, label = TRUE, 
    label.size, large.scale = FALSE, line.width = 1, markers = TRUE, 
    marker.size, max.num, means.plot = FALSE, offset = rep(0.5, 
        4), pch.means = 15, pch.means.size = 1, pos, predictions.mean = NULL, 
    predictions.sample = NULL, ort.lty = 1, specify.bags, specify.classes = 1, 
    strepie = c(1, 1), Title = NULL, Tukey.median = TRUE, Z.means.mat = NULL, 
    straight = TRUE, zoomval = NULL, CLPs = NULL, prediction.regions = NULL) 
{
    if (!is.null(Z.means.mat)) {
        if (ncol(Z.means.mat) == 5) 
            Z.means.mat <- data.frame(Z.means.mat[, 1:5], pch.means.size = pch.means.size, 
                stringsAsFactors = FALSE)
    }
    .orthog.pred.line <- function(x1, y1, x2, y2, val1 = NULL, 
        val2 = NULL, px, py, ort.lty) {
        if (x1 - x2 == 0) 
            gradient <- Inf
        else {
            if (y1 - y2 == 0) 
                gradient <- 0
            else gradient <- (y1 - y2)/(x1 - x2)
        }
        if (gradient == 0) {
            x.star <- px
            y.star <- y1
            prediction <- val1 + (x.star - x1) * (val2 - val1)/(x2 - 
                x1)
        }
        else {
            if (gradient == Inf) {
                x.star <- x1
                y.star <- py
                prediction <- val1 + (y.star - y1) * (val2 - 
                  val1)/(y2 - y1)
            }
            else {
                intcept <- y1 - gradient * x1
                intcept.star <- py + (1/gradient) * px
                x.star <- (intcept.star - intcept)/(gradient + 
                  1/gradient)
                y.star <- intcept + gradient * x.star
                prediction <- val1 + (x.star - x1) * (val2 - 
                  val1)/(x2 - x1)
            }
        }
        lines(x = c(px, x.star), y = c(py, y.star), lty = ort.lty)
        prediction
    }
    .draw.marker <- function(x, y, grad, expand = 1, both.sides = TRUE, 
        col = 1) {
        uin <- par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
        mm <- 1/(uin[1] * 25.4)
        d <- expand * mm
        b <- d * sqrt(1/(1 + grad * grad))
        a <- b * grad
        if (is.infinite(grad)) {
            if (sign(grad) > 0) 
                lines(c(x, x), c(y, y + d), col = col)
            else lines(c(x, x), c(y, y - d), col = col)
        }
        else {
            if (both.sides == FALSE) 
                lines(c(x, x + b), c(y, y + a), col = col)
            else lines(c(x - b, x + b), c(y - a, y + a), col = col)
        }
    }
    .marker.value <- function(x, y, grad, mark, expand = 1, marker.size, 
        col = 1, pos1) {
        uin <- par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
        mm <- 1/(uin[1] * 25.4)
        d <- (expand + 2) * mm
        b <- d * sqrt(1/(1 + grad * grad))
        a <- b * grad
        if (is.infinite(grad)) {
            if (sign(grad) > 0) 
                text(x = x, y = y + d, labels = mark, cex = marker.size, 
                  col = col, pos = 1)
            else text(x = x, y = y - d, labels = mark, cex = marker.size, 
                col = col, pos = 1)
        }
        else {
            text(x = x, y = y, labels = mark, cex = marker.size, 
                col = col, pos = pos1)
        }
    }
    .draw.axis <- function(marker.mat, line.name = NULL, offset = c(0.5, 
        0.5, 0.5, 0.5), pos = "Orthog", axis.col = "black", ax.name.col = "black", 
        ax.name.size = 0.65) {
        if (!is.element(pos[1], c("Hor", "Orthog", "Paral"))) 
            stop("Argument pos must be one of 'Hor','Orthog' or 'Paral' ")
        if (pos[1] == "Hor") {
            par(las = 1)
            adjust <- c(0.5, 1, 0.5, 0)
        }
        if (pos[1] == "Orthog") {
            par(las = 2)
            adjust <- c(0, 0, 0, 0)
        }
        if (pos[1] == "Paral") {
            par(las = 0)
            adjust <- c(0.5, 0.5, 0.5, 0.5)
        }
        marker.mat <- marker.mat[order(marker.mat[, 1]), ]
        x.vals <- marker.mat[, 1]
        y.vals <- marker.mat[, 2]
        marker.vals <- marker.mat[, 3]
        test <- rep(NA, 4)
        test[1] <- min(x.vals) < usr[1]
        test[2] <- max(x.vals) > usr[2]
        test[3] <- min(y.vals) < usr[3]
        test[4] <- max(y.vals) > usr[4]
        if (y.vals[1] == y.vals[length(y.vals)] & x.vals[1] == 
            x.vals[length(x.vals)]) 
            type.line <- "plot.5"
        if (y.vals[1] == y.vals[length(y.vals)] & x.vals[1] != 
            x.vals[length(x.vals)]) {
            gradient <- 0
            intercept <- y.vals[1]
            type.line <- "plot.1"
        }
        if (y.vals[1] != y.vals[length(y.vals)] & x.vals[1] == 
            x.vals[length(x.vals)]) {
            gradient <- Inf
            intercept <- x.vals[1]
            type.line <- "plot.2"
        }
        if (y.vals[1] != y.vals[length(y.vals)] & x.vals[1] != 
            x.vals[length(x.vals)]) {
            gradient <- (y.vals[1] - y.vals[length(y.vals)])/(x.vals[1] - 
                x.vals[length(x.vals)])
            intercept <- y.vals[1] - gradient * x.vals[1]
            y1.ster <- gradient * usr[1] + intercept
            y2.ster <- gradient * usr[2] + intercept
            x1.ster <- (usr[3] - intercept)/gradient
            x2.ster <- (usr[4] - intercept)/gradient
            if (gradient > 0) 
                type.line <- "plot.3"
            if (gradient < 0) 
                type.line <- "plot.4"
        }
        plot.1 <- function(marker.mat, line.name = NULL, offset = NULL, 
            adjust = NULL, axis.col = "black", ax.name.col = "black", 
            ax.name.size = 0.65) {
            abline(h = marker.mat[1, 2], col = axis.col)
            if (!is.null(line.name)) {
                if (ncol(marker.mat) == 3) {
                  if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                    3] < 0) 
                    mtext(text = line.name, side = 4, line = 0 + 
                      offset[4], adj = adjust[4], at = marker.mat[1, 
                      2], col = ax.name.col, cex = ax.name.size)
                  else mtext(text = line.name, side = 2, line = 0 + 
                    offset[2], adj = adjust[2], at = marker.mat[1, 
                    2], col = ax.name.col, cex = ax.name.size)
                }
            }
        }
        plot.2 <- function(marker.mat, line.name = NULL, offset = NULL, 
            adjust = NULL, axis.col = "black", ax.name.col = "black", 
            ax.name.size = 0.65) {
            abline(v = marker.mat[1, 1], col = axis.col)
            if (!is.null(line.name)) {
                if (ncol(marker.mat) == 3) {
                  test <- order(marker.mat[, 2])[c(1, nrow(marker.mat))]
                  if (marker.mat[test[2], 3] - marker.mat[test[1], 
                    3] > 0) 
                    mtext(text = line.name, side = 3, line = 0 + 
                      offset[3], adj = adjust[3], at = marker.mat[1, 
                      1], col = ax.name.col, cex = ax.name.size)
                  else mtext(text = line.name, side = 1, line = 0 + 
                    offset[1], adj = adjust[1], at = marker.mat[1, 
                    1], col = ax.name.col, cex = ax.name.size)
                }
            }
        }
        plot.3 <- function(y1.ster, y2.ster, x1.ster, x2.ster, 
            marker.mat, line.name = NULL, offset = NULL, adjust = NULL, 
            axis.col = "black", ax.name.col = "black", ax.name.size = 0.65) {
            if (y1.ster >= usr[3] & y2.ster >= usr[4]) {
                lines(c(usr[1], x2.ster), c(y1.ster, usr[4]), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 2, line = 0 + 
                        offset[2], adj = adjust[2], at = y1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 3, line = 0 + 
                      offset[3], adj = adjust[3], at = x2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster > usr[3] & y2.ster < usr[4]) {
                lines(c(usr[1], usr[2]), c(y1.ster, y2.ster), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 2, line = 0 + 
                        offset[2], adj = adjust[2], at = y1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 4, line = 0 + 
                      offset[4], adj = adjust[4], at = y2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster < usr[3] & y2.ster > usr[4]) {
                lines(c(x1.ster, x2.ster), c(usr[3], usr[4]), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 1, line = 0 + 
                        offset[1], adj = adjust[1], at = x1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 3, line = 0 + 
                      offset[3], adj = adjust[3], at = x2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster < usr[3] & y2.ster < usr[4]) {
                lines(c(x1.ster, usr[2]), c(usr[3], y2.ster), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 1, line = 0 + 
                        offset[1], adj = adjust[1], at = x1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 4, line = 0 + 
                      offset[4], adj = adjust[4], at = y2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
        }
        plot.4 <- function(y1.ster, y2.ster, x1.ster, x2.ster, 
            marker.mat, line.name = NULL, offset = NULL, adjust = NULL, 
            axis.col = "black", ax.name.col = "black", ax.name.size = 0.65) {
            if (y1.ster > usr[4] & y2.ster > usr[3]) {
                lines(c(x2.ster, usr[2]), c(usr[4], y2.ster), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 3, line = 0 + 
                        offset[3], adj = adjust[3], at = x2.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 4, line = 0 + 
                      offset[4], adj = adjust[4], at = y2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster < usr[4] & y2.ster > usr[3]) {
                lines(c(usr[1], usr[2]), c(y1.ster, y2.ster), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 2, line = 0 + 
                        offset[2], adj = adjust[2], at = y1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 4, line = 0 + 
                      offset[4], adj = adjust[4], at = y2.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster > usr[4] & y2.ster < usr[3]) {
                lines(c(x2.ster, x1.ster), c(usr[4], usr[3]), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) {
                      mtext(text = line.name, side = 3, line = 0 + 
                        offset[3], adj = adjust[3], at = x2.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    }
                    else mtext(text = line.name, side = 1, line = 0 + 
                      offset[1], adj = adjust[1], at = x1.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
            if (y1.ster < usr[4] & y2.ster < usr[3]) {
                lines(c(usr[1], x1.ster), c(y1.ster, usr[3]), 
                  col = axis.col)
                if (!is.null(line.name)) {
                  if (ncol(marker.mat) == 3) {
                    if (marker.mat[1, 3] - marker.mat[nrow(marker.mat), 
                      3] > 0) 
                      mtext(text = line.name, side = 2, line = 0 + 
                        offset[2], adj = adjust[2], at = y1.ster, 
                        col = ax.name.col, cex = ax.name.size)
                    else mtext(text = line.name, side = 1, line = 0 + 
                      offset[1], adj = adjust[1], at = x1.ster, 
                      col = ax.name.col, cex = ax.name.size)
                  }
                }
            }
        }
        plot.5 <- function(marker.mat, axis.col = 1) {
            abline(h = marker.mat[1, 2], col = axis.col)
            abline(v = marker.mat[1, 1], col = axis.col)
        }
        switch(type.line, plot.1 = plot.1(marker.mat = marker.mat, 
            line.name = line.name, offset = offset, adjust = adjust, 
            axis.col = axis.col, ax.name.col = ax.name.col, ax.name.size = ax.name.size), 
            plot.2 = plot.2(marker.mat = marker.mat, line.name = line.name, 
                offset = offset, adjust = adjust, axis.col = axis.col, 
                ax.name.col = ax.name.col, ax.name.size = ax.name.size), 
            plot.3 = plot.3(marker.mat = marker.mat, line.name = line.name, 
                offset = offset, y1.ster = y1.ster, y2.ster = y2.ster, 
                x1.ster = x1.ster, x2.ster = x2.ster, adjust = adjust, 
                axis.col = axis.col, ax.name.col = ax.name.col, 
                ax.name.size = ax.name.size), plot.4 = plot.4(marker.mat = marker.mat, 
                line.name = line.name, offset = offset, y1.ster = y1.ster, 
                y2.ster = y2.ster, x1.ster = x1.ster, x2.ster = x2.ster, 
                adjust = adjust, axis.col = axis.col, ax.name.col = ax.name.col, 
                ax.name.size = ax.name.size), plot.5 = plot.5(marker.mat = marker.mat), 
            axis.col = axis.col)
        list(gradient = gradient, intercept = intercept)
    }
    .axes.plot <- function(ax, z.axes, z.axes.names, predictions.sample, 
        predictions.mean, Z, Z.means.mat, p, ax.col, ax.name.size, 
        markers, label, label.size, marker.size, offset, pos, 
        strepie, ort.lty) {
        if (is.null(ax)) 
            axes <- NULL
        else axes <- (1:p)[ax]
        predictions <- NULL
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
        r.names <- NULL
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
            r.names[i] <- axis.name
            x.invals <- x.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            y.invals <- y.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            tick.labels <- marker.mat[x.vals < usr[2] & x.vals > 
                usr[1] & y.vals < usr[4] & y.vals > usr[3], 3]
            if (length(x.invals) < 2) 
                warning(paste("Less than 2 markers on axis ", 
                  i, ". Increase n.int."))
            uit <- .draw.axis(marker.mat, line.name = axis.name, 
                ax.name.size = ax.name.size, axis.col = ax.col$ax.col[i], 
                ax.name.col = ax.col$ax.name.col[i], offset = offset, 
                pos = pos)
            if (!is.null(predictions.sample)) {
                for (ss in predictions.sample) {
                  predictions[i, paste("s", ss, sep = "")] <- round(.orthog.pred.line(x1 = x.invals[1], 
                    y1 = y.invals[1], x2 = x.invals[length(x.invals)], 
                    y2 = y.invals[length(y.invals)], val1 = tick.labels[1], 
                    val2 = tick.labels[length(x.invals)], px = Z[ss, 
                      1], py = Z[ss, 2], ort.lty = ort.lty), 
                    digits = 4)
                }
            }
            if (!is.null(predictions.mean)) {
                for (ss in predictions.mean) {
                  predictions[i, paste("m", ss, sep = "")] <- round(.orthog.pred.line(x1 = x.invals[1], 
                    y1 = y.invals[1], x2 = x.invals[length(x.invals)], 
                    y2 = y.invals[length(y.invals)], val1 = tick.labels[1], 
                    val2 = tick.labels[length(x.invals)], px = Z.means.mat[ss, 
                      1], py = Z.means.mat[ss, 2], ort.lty = ort.lty), 
                    digits = 4)
                }
            }
            gradient <- uit$gradient
            for (j in 1:length(x.invals)) .draw.marker(x = x.invals[j], 
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
            for (j in 1:length(x.labvals)) .marker.value(x = x.labvals[j], 
                y = y.labvals[j], grad = -1/gradient, mark = tick.labels[j], 
                expand = strepie[2], marker.size = marker.size, 
                col = ax.col$marker.col[i])
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
    .trajectories.plot <- function(z.axes, z.axes.names, p, ax, 
        pos) {
        if (is.null(ax)) 
            axes <- NULL
        else axes <- (1:p)[ax]
        if (!is.element(pos[1], c("Hor", "Orthog", "Paral"))) 
            stop("Argument pos must be one of 'Hor','Orthog' or 'Paral' ")
        if (pos[1] == "Hor") {
            par(las = 1)
            adjust <- c(0.5, 1, 0.5, 0)
        }
        if (pos[1] == "Orthog") {
            par(las = 2)
            adjust <- c(0, 0, 0, 0)
        }
        if (pos[1] == "Paral") {
            par(las = 0)
            adjust <- c(0.5, 0.5, 0.5, 0.5)
        }
        axis.outofbounds <- rep(F, p)
        for (i in axes) {
            if (is.null(z.axes.names)) {
                axis.name <- paste("v", i, sep = "")
            }
            else {
                axis.name <- z.axes.names[i]
            }
            lines(z.axes[[i]][, 1], z.axes[[i]][, 2])
            tick.pos <- (1:nrow(z.axes[[i]]))[z.axes[[i]][, 4] == 
                1]
            tick.pos.voor <- tick.pos - 1
            tick.pos.voor[tick.pos.voor < 1] <- 1
            tick.pos.na <- tick.pos + 1
            tick.pos.na[tick.pos.na > nrow(z.axes[[i]])] <- nrow(z.axes[[i]])
            gradient <- (z.axes[[i]][tick.pos.na, 2] - z.axes[[i]][tick.pos.voor, 
                2])/(z.axes[[i]][tick.pos.na, 1] - z.axes[[i]][tick.pos.voor, 
                1])
            gradient[is.na(gradient)] <- 0
            marker.mat <- z.axes[[i]][tick.pos, 1:3, drop = F]
            x.vals <- marker.mat[, 1]
            y.vals <- marker.mat[, 2]
            x.invals <- x.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            y.invals <- y.vals[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            tick.labels <- marker.mat[x.vals < usr[2] & x.vals > 
                usr[1] & y.vals < usr[4] & y.vals > usr[3], 3]
            gradient <- gradient[x.vals < usr[2] & x.vals > usr[1] & 
                y.vals < usr[4] & y.vals > usr[3]]
            if (length(gradient) == 0) 
                axis.outofbounds[i] <- T
            else {
                for (j in 1:length(x.invals)) {
                  .draw.marker(x.invals[j], y.invals[j], -1/gradient[j])
                  pos1 <- 1
                  if ((-1/gradient[j] > -1) & (-1/gradient[j]) < 
                    1) 
                    pos1 <- 4
                  .marker.value(x.invals[j], y.invals[j], gradient[j], 
                    tick.labels[j], marker.size = marker.size, 
                    pos1 = pos1)
                }
                pos2 <- 2
                if (gradient[length(x.invals)] > -1 & gradient[length(x.invals)] < 
                  1) 
                  pos2 <- 3
                text(x.invals[length(x.invals)], y.invals[length(x.invals)], 
                  axis.name, pos = pos2, cex = ax.name.size)
            }
        }
        write.line <- usr[3] + 0.05 * (usr[4] - usr[3])
        for (i in 1:p) if ((axis.outofbounds[i])) {
            if (is.null(zoomval)) {
                text(usr[1], write.line, paste(paste("Axis", 
                  i), "out of bounds"), pos = 4)
                write.line <- write.line + 0.05 * (usr[4] - usr[3])
            }
        }
    }
    .CLP.plot <- function(CLPs) {
        for (j in 1:length(CLPs)) {
            CLP.mat <- CLPs[[j]]
            text(CLP.mat[, 1], CLP.mat[, 2], dimnames(CLP.mat)[[1]], 
                col = j + 1, cex = 0.7)
        }
    }
    .samples.plot <- function(Z, pch.samples.size, pch.samples, 
        Z.means.mat, specify.classes, label, label.size, class.vec) {
        if (ncol(Z) == 2) 
            Z <- cbind(Z, pch.samples)
        if (ncol(Z) == 3) 
            Z <- cbind(Z, 1)
        if (ncol(Z) == 4) 
            Z <- cbind(Z, 1)
        if (ncol(Z) == 5) 
            Z <- cbind(Z, pch.samples.size)
        if (is.null(dimnames(Z)[[1]])) 
            dimnames(Z) <- list(paste("s", 1:nrow(Z)), NULL)
        classes <- unique(Z[, 3])
        class.vec <- as.numeric(Z[, 3])
        legend.labs <- classes
        for (j in classes) {
            Z.class <- Z[class.vec == j, , drop = FALSE]
            if (label == TRUE) 
                text(Z.class[, 1], Z.class[, 2] - 0.015 * (usr[4] - 
                  usr[3]), labels = dimnames(Z.class)[[1]], cex = label.size)
            Z.class <- data.frame(Z.class[, 1:2], pch.samples.samp = Z.class[, 
                3], col = as.character(Z.class[, 4]), lty = Z.class[1, 
                5], pch.samples.size = Z.class[1, 6], stringsAsFactors = FALSE)
            if (!is.null(specify.classes)) {
                for (i in 1:nrow(Z.class)) points(x = Z.class[i, 
                  1], y = Z.class[i, 2], pch = Z.class[i, 3], 
                  col = as.character(Z.class[i, 4]), cex = Z.class[i, 
                    6])
            }
        }
    }
    .bags.plot <- function(Z, Z.means.mat, specify.bags, class.vec, 
        alpha1, max.num, line.width, Tukey.median, c.hull.n) {
        J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[specify.bags]
        legend.labs.bags <- classes
        for (j in classes) {
            Z.class <- Z[class.vec == j, ]
            chr <- Z.class[1, 3]
            colr <- as.character(Z.class[1, 4])
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
        specify.classes <- NULL
        plot(Z.means.mat[, 1] * exp.factor, Z.means.mat[, 2] * 
            exp.factor, xlim = range(Z.means.mat[, 1] * exp.factor), 
            ylim = range(Z.means.mat[, 2] * exp.factor), xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", type = "n", xaxs = "i", 
            yaxs = "i", asp = 1)
    }
    else {
        if (is.null(zoomval)) 
            plot(Z[, 1] * exp.factor, Z[, 2] * exp.factor, xlim = range(Z[, 
                1] * exp.factor), ylim = range(Z[, 2] * exp.factor), 
                xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
                type = "n", xaxs = "i", yaxs = "i", asp = 1)
        if (is.numeric(zoomval)) 
            plot(Z[, 1] * exp.factor, Z[, 2] * exp.factor, xlim = zoomval[1:2], 
                ylim = zoomval[3:4], xaxt = "n", yaxt = "n", 
                xlab = "", ylab = "", type = "n", xaxs = "i", 
                yaxs = "i", asp = 1)
    }
    usr <- par("usr")
    legend.lab.bags <- NULL
    if (means.plot == TRUE & large.scale == TRUE) {
        if (!is.null(predictions.sample)) {
            predictions.sample <- NULL
            warning("Predictions for sample points only available when large.scale == FALSE \n")
        }
        for (i in 1:nrow(Z.means.mat)) points(x = Z.means.mat[i, 
            1], y = Z.means.mat[i, 2], pch = Z.means.mat[i, 3], 
            col = Z.means.mat[i, 4], cex = Z.means.mat[i, 6])
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
        warning("When means are plotted on a large scale graph no bags or samples are plotted\n")
    }
    if (means.plot == TRUE & large.scale == FALSE) {
        if (!is.null(predictions.mean)) 
            warning("Predictions for class means more suitable when large.scale == TRUE \n")
        for (i in 1:nrow(Z.means.mat)) points(x = Z.means.mat[i, 
            1], y = Z.means.mat[i, 2], pch = Z.means.mat[i, 3], 
            col = Z.means.mat[i, 4], cex = Z.means.mat[i, 6])
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
    }
    if (means.plot == FALSE & large.scale == FALSE) {
        if (!is.null(predictions.mean)) {
            warning("Predictions for class means only available when means.plot == TRUE \n")
            prediction.means <- NULL
        }
        legend.labs <- unique(dimnames(Z.means.mat)[[1]])
    }
    predictions <- NULL
    if (straight) 
        predictions <- .axes.plot(Z = Z, Z.means.mat = Z.means.mat, 
            z.axes = z.axes, z.axes.names = z.axes.names, p = p, 
            ax = ax, ax.name.size = ax.name.size, ax.col = ax.col, 
            markers = markers, label = label, label.size = label.size, 
            marker.size = marker.size, offset = offset, pos = pos, 
            strepie = strepie, predictions.sample = NULL, predictions.mean = predictions.mean, 
            ort.lty = ort.lty)
    else .trajectories.plot(z.axes = z.axes, z.axes.names = z.axes.names, 
        p = p, ax = ax, pos = pos)
    if (!is.null(CLPs)) 
        .CLP.plot(CLPs)
    if (!is.null(specify.classes)) 
        .samples.plot(Z = Z, pch.samples.size = pch.samples.size, 
            pch.samples = pch.samples, Z.means.mat = Z.means.mat, 
            specify.classes = specify.classes, label = label, 
            label.size = label.size, class.vec = class.vec)
    if (!is.null(specify.bags)) 
        legend.lab.bags <- .bags.plot(Z = Z, Z.means.mat = Z.means.mat, 
            specify.bags = specify.bags, class.vec = class.vec, 
            alpha1 = alpha1, max.num = max.num, line.width = line.width, 
            Tukey.median = Tukey.median, c.hull.n = c.hull.n)
    title(main = (ifelse(is.null(Title), paste("Biplot"), Title)))
    if (!is.null(prediction.regions)) {
        for (j in 1:length(prediction.regions)) {
            prediction.mat <- prediction.regions[[j]]
            dev.new()
            if (large.scale) {
                if (nrow(Z.means.mat) == 1) 
                  stop("This option should only be chosen if there are at least two different classes")
                specify.bags <- NULL
                specify.classes <- NULL
                plot(Z.means.mat[, 1] * exp.factor, Z.means.mat[, 
                  2] * exp.factor, xlim = range(Z.means.mat[, 
                  1] * exp.factor), ylim = range(Z.means.mat[, 
                  2] * exp.factor), xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", type = "n", xaxs = "i", yaxs = "i", 
                  asp = 1)
            }
            else {
                if (is.null(zoomval)) 
                  plot(Z[, 1] * exp.factor, Z[, 2] * exp.factor, 
                    xlim = range(Z[, 1] * exp.factor), ylim = range(Z[, 
                      2] * exp.factor), xaxt = "n", yaxt = "n", 
                    xlab = "", ylab = "", type = "n", xaxs = "i", 
                    yaxs = "i", asp = 1)
                if (is.numeric(zoomval)) 
                  plot(Z[, 1] * exp.factor, Z[, 2] * exp.factor, 
                    xlim = zoomval[1:2], ylim = zoomval[3:4], 
                    xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
                    type = "n", xaxs = "i", yaxs = "i", asp = 1)
            }
            for (i in 1:max(prediction.mat[, 3])) points(prediction.mat[prediction.mat[, 
                3] == i, 1], prediction.mat[prediction.mat[, 
                3] == i, 2], col = prediction.mat[prediction.mat[, 
                3] == i, 3] + 5, pch = 16)
        }
    }
    list(legend.lab.bags = legend.lab.bags, predictions = predictions)
}
