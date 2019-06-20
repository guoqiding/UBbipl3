PCAbipl.density.groups <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    X.new.samples = NULL, scaled.mat = TRUE, e.vects = 1:ncol(X), 
    alpha = 0.9, adequacies.print = FALSE, ax = 1:ncol(X), ax.pred.above = NULL, 
    ax.col = list(ax.col = rep(2, ncol(X)), tickmarker.col = rep(3, 
        ncol(X)), marker.col = rep(3, ncol(X))), ax.name.col = rep("black", 
        ncol(X)), ax.name.size = 0.65, ax.type = c("predictive", 
        "interpolative"), bandwidth = NULL, basket.n = 36, basket.beta = 0.4, 
    between = c(1, -1, 0, 1), between.columns = -1, char.legend.size = c(1.2, 
        0.7), c.hull.n = 10, colours = c(4:12, 3:1), colour.scheme = NULL, 
    colours.means = NULL, colours.density = matrix(rep(c("white", 
        "red"), ncol(G)), byrow = T, nrow = ncol(G)), columns = 1, 
    correlation.biplot = FALSE, cuts.density = 50, dens.ax.cex = 0.6, 
    dens.ax.tcl = -0.2, dens.mgp = c(0, -0.25, 0), draw.densitycontours = FALSE, 
    exp.factor = 1.2, label = TRUE, label.size = 0.6, large.scale = FALSE, 
    layout.heights = c(100, 10), legend.type = c(means = TRUE, 
        samples = TRUE, bags = FALSE), line.length = c(1, 1), 
    line.size = 2.5, line.type = 1:ncol(G), line.width = 1, markers = TRUE, 
    marker.size = 0.5, max.num = 2500, means.plot = TRUE, n.int = rep(5, 
        ncol(X)), offset = rep(0.5, 4), offset.m = rep(0.001, 
        ncol(X)), parlegendmar = c(0, 6, 2.5, 6), parlegendmar.dens = c(2, 
        5, 0, 5), parplotmar = rep(3, 4), pch.means = 0:10, pch.means.size = 1, 
    pch.samples = 0:10, pch.samples.size = 1, pos = c("Orthog", 
        "Hor", "Paral"), pos.m = rep(1, ncol(X)), predictivity.print = FALSE, 
    quality.print = FALSE, select.origin = FALSE, side.label = rep("right", 
        ncol(X)), smooth.n = 100, specify.bags = dimnames(G)[[2]], 
    specify.classes = dimnames(G)[[2]], specify.density.classes = dimnames(G)[[2]], 
    specify.beta.baskets = NULL, Title = NULL, Tukey.median = TRUE) 
{
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    ax.type <- ax.type[1]
    pos <- pos[1]
    if (is.null(G)) {
        G <- matrix(indmat(rep(1, nrow(X))), ncol = 1)
        dimnames(G) <- list(1:nrow(X), "AllData")
    }
    old.par <- par(no.readonly = TRUE)
    if (adequacies.print & predictivity.print) 
        stop("adequacies.print and predictivity.print cannot both be set to True")
    if (!is.null(G)) 
        if ((length(specify.bags) > ncol(G)) | (length(specify.classes) > 
            ncol(G)) | (length(specify.density.classes) > ncol(G))) 
            return(cat("Number of specified bags or specified classes or specified density classes must not be larger than the number of different classes\n"))
    layout(matrix(1:2, ncol = 1), heights = layout.heights)
    par(pty = "s", mar = parplotmar, cex = 1)
    on.exit(par(old.par))
    X <- as.matrix(X)
    unscaled.X <- X
    means <- apply(X, 2, mean)
    sd <- sqrt(apply(X, 2, var))
    if (scaled.mat) 
        X <- scale(X)
    else {
        X <- scale(X, scale = FALSE)
        sd <- rep(1, ncol(X))
    }
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    n.groups <- apply(G, 2, sum)
    if (J > length(pch.means)) 
        stop(paste("Increase size of pch.means argument"))
    if (J > length(pch.samples)) 
        stop(paste("Increase size of pch.samples argument"))
    if (J > length(colours)) 
        stop(paste("Increase size of colours argument"))
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste(1:n), paste("V", 1:p, sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste(1:n)
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
    if (!is.null(G)) {
        if (nrow(G) != n) 
            stop("number of rows of X and G differ")
        if (is.null(dimnames(G))) 
            dimnames(G) <- list(NULL, paste("class", 1:J, sep = ""))
        if (length(dimnames(G)[[2]]) == 0) 
            dimnames(G)[[2]] <- paste("class", 1:J, sep = "")
        if (is.numeric(specify.bags)) 
            specify.bags <- dimnames(G)[[2]][specify.bags]
        if (!is.null(specify.bags)) {
            if (!is.numeric(specify.bags)) {
                test <- match(specify.bags, unique(as.character(dimnames(G)[[2]])), 
                  nomatch = 0)
                if (sum(test != 0) < length(specify.bags)) 
                  stop(" One or more of specified classes non-existent.  Check spelling")
            }
        }
        if (!is.null(specify.bags)) {
            select.numeric.bags <- rep(0, length(specify.bags))
            for (k in 1:length(specify.bags)) select.numeric.bags[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                specify.bags[k]]
        }
        else {
            if (is.null(specify.bags)) 
                select.numeric.bags <- NULL
        }
        if (is.numeric(specify.beta.baskets)) 
            specify.beta.baskets <- dimnames(G)[[2]][specify.beta.baskets]
        if (!is.null(specify.beta.baskets)) {
            if (!is.numeric(specify.beta.baskets)) {
                test <- match(specify.beta.baskets, unique(as.character(dimnames(G)[[2]])), 
                  nomatch = 0)
                if (sum(test != 0) < length(specify.beta.baskets)) 
                  stop(" One or more of specified classes non-existent.  Check spelling")
            }
        }
        if (!is.null(specify.beta.baskets)) {
            select.numeric.beta.baskets <- rep(0, length(specify.beta.baskets))
            for (k in 1:length(specify.beta.baskets)) select.numeric.beta.baskets[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                specify.beta.baskets[k]]
        }
        else {
            if (is.null(specify.beta.baskets)) 
                select.numeric.beta.baskets <- NULL
        }
        if (is.numeric(specify.classes)) 
            specify.classes <- dimnames(G)[[2]][specify.classes]
        if (!is.null(specify.classes)) {
            if (!is.numeric(specify.classes)) {
                test <- match(specify.classes, unique(as.character(dimnames(G)[[2]])), 
                  nomatch = 0)
                if (sum(test != 0) < length(specify.classes)) 
                  stop(" One or more of specified classes non-existent.  Check spelling")
            }
        }
        if (!is.null(specify.classes)) {
            select.numeric.classes <- rep(0, length(specify.classes))
            for (k in 1:length(specify.classes)) select.numeric.classes[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                specify.classes[k]]
        }
        else {
            if (is.null(specify.classes)) 
                select.numeric.classes <- NULL
        }
        if (is.numeric(specify.density.classes)) 
            specify.density.classes <- dimnames(G)[[2]][specify.density.classes]
        if (!is.null(specify.density.classes)) {
            if (!is.numeric(specify.density.classes)) {
                test <- match(specify.density.classes, unique(as.character(dimnames(G)[[2]])), 
                  nomatch = 0)
                if (sum(test != 0) < length(specify.density.classes)) 
                  stop(" One or more of specified density classes non-existent.  Check spelling")
            }
        }
        if (!is.null(specify.density.classes)) {
            select.numeric.density.classes <- rep(0, length(specify.density.classes))
            for (k in 1:length(specify.density.classes)) select.numeric.density.classes[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                specify.density.classes[k]]
        }
        else {
            if (is.null(specify.density.classes)) 
                select.numeric.density.classes <- NULL
        }
        dimnames(G)[[2]] <- paste(dimnames(G)[[2]], "; n = ", 
            n.groups, sep = "")
        if (length(dimnames(G)[[2]]) == 1) 
            class.vec <- rep(dimnames(G)[[2]], nrow(X))
        else {
            class.vec <- apply(t(apply(G, 1, function(x) x == 
                max(x))), 1, function(s, G) dimnames(G)[[2]][s], 
                G = G)
        }
    }
    svd.out <- svd(t(X) %*% X)
    V.mat <- svd.out$u
    Vr <- svd.out$u[, e.vects[1:2]]
    eigval <- svd.out$d
    lambda.mat <- diag(eigval)
    eigval.r <- eigval[e.vects[1:2]]
    lambda.r.mat <- diag(eigval.r)
    fit.predictivity.mat <- diag(diag(Vr %*% lambda.r.mat %*% 
        t(Vr))) %*% solve(diag(diag(V.mat %*% lambda.mat %*% 
        t(V.mat))))
    fit.predictivity <- round(diag(fit.predictivity.mat), digits = 3)
    names(fit.predictivity) <- dimnames(X)[[2]]
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
        "%")
    fit.adequacy <- round(diag(Vr %*% t(Vr)), digits = 3)
    names(fit.adequacy) <- dimnames(X)[[2]]
    if (adequacies.print) 
        dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", fit.adequacy, 
            ")", sep = "")
    if (predictivity.print) 
        dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", round(fit.predictivity, 
            digits = 2), ")", sep = "")
    if (!is.null(G)) {
        class.means.mat <- solve(t(G) %*% G) %*% t(G) %*% unscaled.X
    }
    if (!is.null(ax.pred.above)) 
        ax <- (1:p)[fit.predictivity >= ax.pred.above]
    if (correlation.biplot) {
        lambda.r <- diag(svd(t(X) %*% X)$d[1:2])
        Z <- sqrt(n - 1) * X %*% Vr %*% (sqrt(solve(lambda.r)))
        if (!is.null(class.means.mat)) {
            Z.means.mat <- sqrt(n - 1) * scale(class.means.mat, 
                means, sd) %*% Vr %*% (sqrt(solve(lambda.r)))
            Z.means.mat <- data.frame(Z.means.mat, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                cex.means = pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
        }
    }
    else {
        Z <- X %*% Vr
        if (!is.null(G)) {
            Z.means.mat <- scale(class.means.mat, means, sd) %*% 
                Vr
            Z.means.mat <- data.frame(Z.means.mat, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                cex.means = pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
        }
    }
    Z <- data.frame(Z, pch.samples = pch.samples[1], colr = as.character(colours[1]), 
        line.type = line.type[1], stringsAsFactors = FALSE)
    if (J > 0) {
        for (j in 1:J) {
            Z[G[, j] == 1, 4] <- as.character(colours[j])
            Z[G[, j] == 1, 3] <- pch.samples[j]
            Z[G[, j] == 1, 5] <- line.type[j]
        }
        if (is.null(specify.bags)) {
            line.type <- NULL
            Z <- Z[, 1:4]
        }
    }
    if (ax.type == "predictive") 
        axes.rows <- 1/(diag(Vr %*% t(Vr))) * Vr
    else {
        if (ax.type == "interpolative") 
            axes.rows <- Vr
        else stop("ax.type is either 'interpolative' or 'predictive' (the default)")
    }
    if (correlation.biplot) {
        axes.rows <- (sqrt(n - 1)/(diag(Vr %*% lambda.r %*% t(Vr)))) * 
            Vr %*% sqrt(lambda.r)
    }
    z.axes <- lapply(1:p, function(j, unscaled.X, means, sd, 
        axes.rows, n.int) {
        number.points <- 100
        std.markers <- pretty(unscaled.X[, j], n = n.int[j])
        std.range <- c(min(std.markers), max(std.markers))
        std.markers.min <- std.markers - (std.range[2] - std.range[1])
        std.markers.max <- std.markers + (std.range[2] - std.range[1])
        std.markers <- c(std.markers, std.markers.min, std.markers.max)
        interval <- (std.markers - means[j])/sd[j]
        axis.vals <- seq(from = min(interval), to = max(interval), 
            length = number.points)
        axis.vals <- sort(unique(c(axis.vals, interval)))
        number.points <- length(axis.vals)
        axis.points <- matrix(0, nrow = number.points, ncol = 4)
        axis.points[, 1] <- axis.vals * axes.rows[j, 1]
        axis.points[, 2] <- axis.vals * axes.rows[j, 2]
        axis.points[, 3] <- axis.vals * sd[j] + means[j]
        axis.points[, 4] <- 0
        for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
            3] - std.markers) == 0)) 
            axis.points[i, 4] <- 1
        return(axis.points)
    }, unscaled.X = unscaled.X, means = means, sd = sd, axes.rows = axes.rows, 
        n.int = n.int)
    uit <- drawbipl.bagalpha.density.groups(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
        p = p, alpha = alpha, basket.beta = basket.beta, basket.n = basket.n, 
        ax = ax, ax.col = ax.col, ax.name.col = ax.name.col, 
        ax.name.size = ax.name.size, bandwidth = bandwidth, c.hull.n = c.hull.n, 
        class.vec = class.vec, colours.density = colours.density, 
        cuts.density = cuts.density, dens.ax.cex = dens.ax.cex, 
        dens.ax.tcl = dens.ax.tcl, dens.mgp = dens.mgp, draw.densitycontours = draw.densitycontours, 
        exp.factor = exp.factor, label = label, label.size = label.size, 
        line.width = line.width, large.scale = large.scale, layout.heights = layout.heights, 
        marker.size = marker.size, markers = markers, max.num = max.num, 
        means.plot = means.plot, offset = offset, offset.m = offset.m, 
        parlegendmar.dens = parlegendmar.dens, pch.means = pch.means, 
        pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
        pos = pos, pos.m = pos.m, smooth.n = smooth.n, specify.bags = select.numeric.bags, 
        specify.classes = select.numeric.classes, specify.density.classes = select.numeric.density.classes, 
        specify.beta.baskets = select.numeric.beta.baskets, side.label = side.label, 
        strepie = line.length, Title = Title, Tukey.median = Tukey.median, 
        Z.means.mat = Z.means.mat)
    if (select.origin) {
        if (ax.type != "predictive") 
            stop("Select.origin only available for predictive biplot axes \n")
        cat("Move cursor where you want origin and press left mouse button \n")
        flush.console()
        origin.pos <- locator(1)
    }
    par(pty = "m", mar = parlegendmar.dens)
    draw.density.rect(levels.rect = uit$levels.rect, col.use = uit$col.use, 
        dens.ax.cex = dens.ax.cex, dens.mgp = dens.mgp, dens.ax.tcl = dens.ax.tcl)
    if (!(is.null(X.new.samples))) {
        pch.new <- paste("N", 1:nrow(X.new.samples), sep = "")
        if (correlation.biplot) {
            Z.new <- sqrt(n - 1) * scale(X.new.samples, means, 
                sd) %*% Vr %*% (sqrt(solve(lambda.r)))
        }
        else Z.new <- scale(X.new.samples, means, sd) %*% Vr
        text(Z.new, labels = pch.new, cex = pch.samples.size, 
            col = "black")
    }
    if (select.origin) {
        z.axes <- lapply(1:p, function(j, unscaled.X, means, 
            sd, axes.rows, n.int, orthog.transx, orthog.transy, 
            p) {
            phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                axes.rows %*% c(orthog.transx[j], orthog.transy[j])
            number.points <- 100
            std.markers <- pretty(unscaled.X[, j], n = n.int[j])
            std.range <- c(min(std.markers), max(std.markers))
            std.markers.min <- std.markers - (std.range[2] - 
                std.range[1])
            std.markers.max <- std.markers + (std.range[2] - 
                std.range[1])
            std.markers <- c(std.markers, std.markers.min, std.markers.max)
            interval <- (std.markers - means[j])/sd[j]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            number.points <- length(axis.vals)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 1]
            axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 2]
            axis.points[, 3] <- axis.vals * sd[j] + means[j]
            axis.points[, 4] <- 0
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, means = means, sd = sd, axes.rows = axes.rows, 
            n.int = n.int, orthog.transx = rep(origin.pos$x, 
                p), orthog.transy = rep(origin.pos$y, p), p = p)
        dev.new()
        layout(matrix(1:2, ncol = 1), heights = layout.heights)
        par(pty = "s", mar = parplotmar, cex = 1)
        drawbipl.bagalpha.density(Z = Z, z.axes = z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.col = ax.name.col, 
            ax.name.size = ax.name.size, alpha = alpha, basket.beta = basket.beta, 
            basket.n = basket.n, ellipse.kappa = ellipse.kappa, 
            ellipse.alpha = ellipse.alpha, pch.means = pch.means, 
            pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
            specify.bags = select.numeric.bags, specify.ellipses = select.numeric.ellipses, 
            label = label, markers = markers, Title = Title, 
            means.plot = means.plot, large.scale = large.scale, 
            specify.classes = select.numeric.classes, specify.beta.baskets = select.numeric.beta.baskets, 
            Tukey.median = Tukey.median, Z.means.mat = Z.means.mat, 
            offset = offset, offset.m = offset.m, pos = pos, 
            pos.m = pos.m, strepie = line.length, max.num = max.num, 
            c.hull.n = c.hull.n, marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, side.label = side.label, smooth.n = smooth.n, 
            cuts.density = cuts.density, colours.density = colours.density, 
            draw.densitycontours = draw.densitycontours, parlegendmar.dens = parlegendmar.dens, 
            dens.ax.cex = dens.ax.cex, dens.ax.tcl = dens.ax.tcl, 
            dens.mgp = dens.mgp, bandwidth = bandwidth)
        par(pty = "m", mar = parlegendmar.dens)
        draw.density.rect(levels.rect = uit$levels.rect, col.use = uit$col.use, 
            dens.ax.cex = dens.ax.cex, dens.mgp = dens.mgp, dens.ax.tcl = dens.ax.tcl)
        if (!(is.null(X.new.samples))) {
            if (is.null(pch.new)) {
                pch.new <- paste("N", 1:nrow(X.new.samples), 
                  sep = "")
                if (correlation.biplot) {
                  Z.new <- sqrt(n - 1) * scale(X.new.samples, 
                    means, sd) %*% Vr %*% (sqrt(solve(lambda.r)))
                }
                else Z.new <- scale(X.newsamples, means, sd) %*% 
                  Vr
                text(Z.new, labels = pch.new, cex = pch.samples.size, 
                  col = "black")
            }
            else {
                if (correlation.biplot) {
                  Z.new <- sqrt(n - 1) * scale(X.new.samples, 
                    means, sd) %*% Vr %*% (sqrt(solve(lambda.r)))
                }
                else Z.new <- scale(X.new.samples, means, sd) %*% 
                  Vr
                points(Z.new, pch = pch.new, cex = pch.samples.size, 
                  col = pch.new.cols)
            }
        }
    }
    if (!is.null(uit$legend.lab.bags)) 
        dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
    if (any(legend.type)) {
        dev.new()
        blegend.colchar(quality.print = quality.print, QualityOfDisplay = fit.quality, 
            classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J]), 
            pch.samples = as.vector(pch.samples[1:J]), colours = as.vector(colours[1:J]), 
            colours.means = colours.means, line.type = as.vector(line.type[1:J]), 
            line.width = line.width, pch.samples.size = char.legend.size, 
            legend.type = legend.type, between = between, columns = columns, 
            parlegendmar = parlegendmar, line.size = line.size, 
            between.columns = between.columns)
    }
    list(Z = Z, Eigenvectors = Vr, e.vals = eigval, quality = fit.quality, 
        adequacy = fit.adequacy, predictivity = fit.predictivity)
}
