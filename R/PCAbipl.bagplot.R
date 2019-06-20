PCAbipl.bagplot <-
function (X = Ocotea.data[, 3:8], G = NULL, X.new = NULL, scaled.mat = TRUE, 
    e.vects = 1:ncol(X), dim.biplot = 2, adequacies.print = FALSE, 
    ax = 1:ncol(X), ax.name.col = rep("black", ncol(X)), ax.name.size = 0.65, 
    ax.type = c("predictive", "interpolative"), ax.col = list(ax.col = rep("grey", 
        ncol(X)), tickmarker.col = rep("grey", ncol(X)), marker.col = rep("black", 
        ncol(X))), between = c(1, -1, 0, 1), between.columns = -1, 
    char.legend.size = c(1.2, 0.7), c.hull.n = 10, colour.scheme = NULL, 
    colours = c(4:12, 3:1), columns = 1, constant = 0.1, correlation.biplot = FALSE, 
    exp.factor = 1.2, label = TRUE, label.size = 0.6, legend.type = c(means = TRUE, 
        samples = TRUE, bags = FALSE), line.length = c(1, 1), 
    line.size = 2.5, line.type = 1:ncol(G), line.width = 1, markers = TRUE, 
    marker.size = 0.5, max.num = 2500, means.plot = TRUE, means.only = FALSE, 
    n.int = rep(5, ncol(X)), offset = rep(0.5, 4), ort.lty = 1, 
    parlegendmar = c(3, 1, 3, 1), parplotmar = rep(3, 4), pch.means = 0:10, 
    pch.means.size = 1, pch.samples = 0:10, pch.samples.size = 1, 
    col.mean = "black", pos = c("Orthog", "Hor", "Paral"), predictions.mean = NULL, 
    predictions.sample = NULL, quality.print = FALSE, resolutions.print = FALSE, 
    specify.bags = dimnames(G)[[2]], specify.classes = dimnames(G)[[2]], 
    text.width.mult = 1, Title = NULL, Tukey.median = TRUE, colr.1 = "grey", 
    colr.3 = "lightgrey", ...) 
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
    if (!(dim.biplot == 2)) 
        stop("Bagplots are drawn only when dim.biplot = 2 \n")
    old.par <- par(no.readonly = TRUE)
    if (!is.null(predictions.mean) & !is.null(predictions.sample)) 
        stop("Either predictions.mean or predictions.sample or both must be NULL \n")
    if (adequacies.print & resolutions.print) 
        stop("adequacies.print and resolutions.print cannot both be set to True")
    if (!is.null(G)) 
        if ((length(specify.bags) > ncol(G)) | (length(specify.classes) > 
            ncol(G))) 
            return(cat("Number of specified bags or specified classes must not be larger than the number of different classes\n"))
    if (!is.null(colour.scheme)) 
        palette(colour.scheme)
    par(pty = "s", mar = parplotmar)
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
            stop("Class names must be supplied, not numeric values")
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
        if (is.numeric(specify.classes)) 
            stop("Class names must be supplied, not numeric values")
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
    fit.resolution.mat <- diag(diag(Vr %*% lambda.r.mat %*% t(Vr))) %*% 
        solve(diag(diag(V.mat %*% lambda.mat %*% t(V.mat))))
    fit.resolution <- round(diag(fit.resolution.mat), digits = 3)
    names(fit.resolution) <- dimnames(X)[[2]]
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
        "%")
    fit.adequacy <- round(diag(Vr %*% t(Vr)), digits = 3)
    names(fit.adequacy) <- dimnames(X)[[2]]
    if (adequacies.print) 
        dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", fit.adequacy, 
            ")", sep = "")
    if (resolutions.print) 
        dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", round(fit.resolution, 
            digits = 2), ")", sep = "")
    if (!is.null(G)) {
        class.means.mat <- as.matrix(solve(t(G) %*% G) %*% t(G) %*% 
            unscaled.X, ncol = ncol(unscaled.X))
    }
    if (correlation.biplot) {
        lambda.r <- diag(svd(t(X) %*% X)$d[1:2])
        Z <- sqrt(n - 1) * X %*% Vr %*% (sqrt(solve(lambda.r)))
        if (!is.null(class.means.mat)) {
            Z.means.mat <- sqrt(n - 1) * scale(class.means.mat, 
                means, sd) %*% Vr %*% (sqrt(solve(lambda.r)))
            Z.means.mat <- data.frame(Z.means.mat, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                pch.means.size, stringsAsFactors = FALSE)
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
                pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
        }
    }
    Z <- data.frame(Z, pch.sampl = pch.samples[1], colr = as.character(colours[1]), 
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
    {
        graphics.off()
        J <- nrow(Z.means.mat)
        classes <- unique(dimnames(Z.means.mat)[[1]])[select.numeric.classes]
        legend.labs.bags <- classes
        for (j in classes) {
            dev.new()
            uit <- drawbipl.bagalpha(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
                p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
                pch.means = pch.means, pch.means.size = 0, pch.samples.size = 0, 
                specify.bags = NULL, label = label, markers = markers, 
                Title = "", means.plot = TRUE, large.scale = FALSE, 
                specify.classes = NULL, Tukey.median = FALSE, 
                c.hull.n = c.hull.n, Z.means.mat = Z.means.mat, 
                offset = offset, pos = pos, strepie = line.length, 
                max.num = max.num, marker.size = marker.size, 
                label.size = label.size, exp.factor = exp.factor, 
                line.width = line.width, class.vec = class.vec, 
                predictions.sample = predictions.sample, predictions.mean = predictions.mean, 
                ort.lty = ort.lty)
            Z.class <- data.frame(Z[class.vec == j, ], stringsAsFactors = FALSE)
            chr <- Z.class[1, 3]
            colr <- Z.class[1, 4]
            Z.class.mean <- Z.means.mat[match(j, classes), , 
                drop = FALSE]
            flush.console()
            cat(paste("bag for class ", j, " with ", nrow(Z.class), 
                " samples in class ", j, sep = ""), "\n")
            my.bagplot(Z.class[, 1], Z.class[, 2], colr.1 = colr.1, 
                colr.3 = colr.3, colr.2 = colr, max.num = max.num, 
                ...)
            if (means.plot) {
                points(x = Z.class.mean[1, 1], y = Z.class.mean[1, 
                  2], pch = Z.class.mean[1, 3], col = col.mean, 
                  cex = pch.means.size)
            }
            if (means.only == FALSE) {
                usr <- par("usr")
                if (is.null(dimnames(Z.class))) 
                  dimnames(Z.class) <- list(paste("s", 1:nrow(Z.class), 
                    sep = ""), NULL)
                if (length(dimnames(Z.class)[[1]]) == 0) 
                  dimnames(Z.class) <- list(paste("s", 1:nrow(Z.class), 
                    sep = ""), dimnames(Z.class)[[2]])
                if (label == TRUE) 
                  text(Z.class[, 1], Z.class[, 2] - 0.015 * (usr[4] - 
                    usr[3]), labels = dimnames(Z.class)[[1]], 
                    cex = 0.65)
                points(Z.class[, 1], Z.class[, 2], col = as.vector(colr), 
                  cex = pch.samples.size)
            }
        }
    }
    if (!(is.null(X.new))) {
        pch.new <- paste("N", 1:nrow(X.new), sep = "")
        if (correlation.biplot) {
            Z.new <- sqrt(n - 1) * scale(X.new, means, sd) %*% 
                Vr %*% (sqrt(solve(lambda.r)))
        }
        else Z.new <- scale(X.new, means, sd) %*% Vr
        text(Z.new, labels = pch.new, cex = pch.samples.size, 
            col = "black")
    }
    if (!is.null(uit$legend.lab.bags)) 
        dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
    if (any(legend.type)) {
        dev.new()
        blegend.colchar(quality.print = quality.print, QualityOfDisplay = fit.quality, 
            classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J]), 
            pch.samples = as.vector(pch.samples[1:J]), colours = as.vector(colours[1:J]), 
            line.type = as.vector(line.type[1:J]), line.width = line.width, 
            pch.samples.size = char.legend.size, legend.type = legend.type, 
            between = between, columns = columns, parlegendmar = parlegendmar, 
            line.size = line.size, between.columns = between.columns, 
            text.width.mult = text.width.mult)
    }
    list(Z = Z, Z.axes = z.axes, V = V.mat, Eigenvectors = Vr, 
        e.vals = eigval, quality = fit.quality, adequacy = fit.adequacy, 
        resolution = fit.resolution, predictions = uit$predictions)
}
