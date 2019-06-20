PCAbipl.OneDim <-
function (X = Ocotea.data[, 3:8], G = indmat(rep(1, nrow(X))), 
    X.new.samples = NULL, scaled.mat = FALSE, e.vects = 1:ncol(X), 
    adequacies.print = FALSE, alpha = 0.95, alpha.3d = 0.7, aspect.3d = "iso", 
    ax.col.3d = "black", ax = 1:ncol(X), ax.name.col = rep("black", 
        ncol(X)), ax.name.size = 0.65, ax.type = c("predictive", 
        "interpolative"), ax.col = list(ax.col = rep(8, ncol(X)), 
        tickmarker.col = rep(8, ncol(X)), marker.col = rep(1, 
            ncol(X))), between = c(1, -1, 0, 1), between.columns = -1, 
    cex.3d = 0.6, char.legend.size = c(1.2, 0.7), c.hull.n = 10, 
    colours = c(4:12, 3:1), colour.scheme = NULL, col.plane.3d = "lightgrey", 
    col.text.3d = "black", columns = 1, constant = 0.1, correlation.biplot = FALSE, 
    density.plot = FALSE, exp.factor = 1.2, factor.x = 2, factor.y = 2, 
    font.3d = 2, ID.labs = FALSE, ID.3d = 1:nrow(X), label = TRUE, 
    label.size = 0.6, large.scale = FALSE, legend.type = c(means = FALSE, 
        samples = TRUE, bags = FALSE), line.length = c(1, 1), 
    line.size = 2.5, line.type = 1:ncol(G), line.width = 1, markers = TRUE, 
    marker.size = 0.5, max.num = 2500, means.plot = TRUE, n.int = rep(5, 
        ncol(X)), offset = rep(0.5, 4), offset.m = rep(0, ncol(X)), 
    ort.lty = 1, parlegendmar = c(3, 1, 3, 1), parplotmar = rep(3, 
        4), pch.means = 0:10, pch.means.size = 1, pch.samples = 0:10, 
    pch.samples.size = 1, pos = c("Orthog", "Hor", "Paral"), 
    pos.m = rep(1, ncol(X)), side.label = rep("right", ncol(X)), 
    predictions.mean = NULL, predictions.sample = NULL, predictivity.print = FALSE, 
    quality.print = FALSE, reflect = diag(1), size.ax.3d = 0.5, 
    size.means.3d = 10, size.points.3d = 5, specify.bags = NULL, 
    specify.classes = dimnames(G)[[2]], text.width.mult = 1, 
    Title = NULL, Titles.3d = c("", "", "x", "y", "z"), Tukey.median = TRUE) 
{
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    ax.type <- ax.type[1]
    pos <- pos[1]
    dim.biplot <- 1
    if (!is.null(specify.bags)) 
        stop("Set specify.bags equal to NULL. Alpha bags are not constructed for 1D biplots. \n")
    old.par <- par(no.readonly = TRUE)
    if (adequacies.print & predictivity.print) 
        stop("adequacies.print and predictivity.print cannot both be set to True")
    if (!is.null(G)) 
        if ((length(specify.bags) > ncol(G)) | (length(specify.classes) > 
            ncol(G))) 
            return(cat("Number of specified bags or specified classes must not be larger than the number of different classes\n"))
    par(pty = "m", mar = parplotmar)
    on.exit(par(old.par))
    reflect.mat <- diag(dim.biplot)
    if ((reflect == "x") || (reflect == "y")) 
        if (dim.biplot == 1 && reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    rotate.mat <- diag(dim.biplot)
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
    V.1 <- svd.out$u[, e.vects[1], drop = FALSE] %*% rotate.mat %*% 
        reflect.mat
    eigval <- svd.out$d
    lambda.mat <- diag(eigval)
    eigval.1 <- eigval[e.vects[1]]
    lambda.1.mat <- as.matrix(eigval.1)
    fit.predictivity.mat <- diag(diag(V.1 %*% lambda.1.mat %*% 
        t(V.1))) %*% solve(diag(diag(V.mat %*% lambda.mat %*% 
        t(V.mat))))
    fit.predictivity <- round(diag(fit.predictivity.mat), digits = 3)
    names(fit.predictivity) <- dimnames(X)[[2]]
    fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]])/sum(eigval)) * 
        100, digits = 2), "%")
    fit.adequacy <- round(diag(V.1 %*% t(V.1)), digits = 3)
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
    if (correlation.biplot) {
        lambda.1 <- diag(svd(t(X) %*% X)$d[1])
        Z <- sqrt(n - 1) * X %*% V.1 %*% (sqrt(solve(lambda.1)))
        if (!is.null(class.means.mat)) {
            Z.means.mat <- sqrt(n - 1) * scale(class.means.mat, 
                means, sd) %*% V.1 %*% (sqrt(solve(lambda.1)))
            Z.means.mat <- data.frame(Z.means.mat, 0, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                pch.means.size = pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
        }
    }
    else {
        Z <- X %*% V.1
        if (!is.null(G)) {
            Z.means.mat <- scale(class.means.mat, means, sd) %*% 
                V.1
            Z.means.mat <- data.frame(Z.means.mat, 0, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                pch.means.size = pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
        }
    }
    Z <- data.frame(Z, 0, pch.sample = pch.samples[1], colr = as.character(colours[1]), 
        line.type = line.type[1], cex.sampl = pch.samples.size, 
        stringsAsFactors = FALSE)
    if (J > 0) {
        for (j in 1:J) {
            Z[G[, j] == 1, 4] <- colours[j]
            Z[G[, j] == 1, 3] <- pch.samples[j]
            Z[G[, j] == 1, 5] <- line.type[j]
        }
        if (is.null(specify.bags)) {
            line.type <- NULL
            Z <- Z[, c(1:4, 6)]
        }
    }
    if (ax.type == "predictive") 
        axes.rows <- 1/(diag(V.1 %*% t(V.1))) * V.1
    else {
        if (ax.type == "interpolative") 
            axes.rows <- V.1
        else stop("ax.type is either 'interpolative' or 'predictive' (the default)")
    }
    if (correlation.biplot) {
        axes.rows <- (sqrt(n - 1)/(diag(V.1 %*% lambda.1 %*% 
            t(V.1)))) * V.1 %*% sqrt(lambda.1)
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
        axis.points[, 2] <- 0
        axis.points[, 3] <- axis.vals * sd[j] + means[j]
        axis.points[, 4] <- 0
        for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
            3] - std.markers) == 0)) 
            axis.points[i, 4] <- 1
        return(axis.points)
    }, unscaled.X = unscaled.X, means = means, sd = sd, axes.rows = axes.rows, 
        n.int = n.int)
    uit <- drawbipl.onedim.bagalpha(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
        p = p, ax = ax, ax.col = ax.col, ax.name.col = ax.name.col, 
        ax.name.size = ax.name.size, alpha = alpha, pch.means = pch.means, 
        pch.means.size = pch.means.size, pch.samples.size = pch.samples.size, 
        specify.bags = select.numeric.bags, label = label, markers = markers, 
        Title = Title, means.plot = means.plot, large.scale = large.scale, 
        specify.classes = select.numeric.classes, Tukey.median = Tukey.median, 
        Z.means.mat = Z.means.mat, offset = offset, offset.m = offset.m, 
        ort.lty = ort.lty, pos = pos, pos.m = pos.m, side.label = side.label, 
        strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
        marker.size = marker.size, label.size = label.size, exp.factor = exp.factor, 
        line.width = line.width, class.vec = class.vec, predictions.sample = predictions.sample, 
        predictions.mean = predictions.mean, constant = constant, 
        density.plot = density.plot)
    if (!(is.null(X.new.samples))) {
        Z.new <- scale(X.new.samples, X.means, scale = FALSE) %*% 
            Br %*% rotate.mat %*% reflect.mat
        if (is.null(pch.new)) {
            pch.new <- paste("N", 1:nrow(X.new.samples), sep = "")
            points(Z.new, pch = pch.new, cex = pch.samples.size, 
                col = pch.new.cols)
            text(Z.new, labels = pch.new, cex = pch.samples.size, 
                col = "black")
        }
        else {
            points(Z.new, pch = pch.new, cex = pch.samples.size, 
                col = pch.new.cols)
            text(Z.new, labels = pch.new.labels, pos = 1, cex = 0.75, 
                col = pch.new.cols)
        }
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
            line.size = line.size, between.columns = between.columns)
    }
    list(Z = Z, Z.axes = z.axes, V = V.mat, Eigenvectors = V.1, 
        e.vals = eigval, PCA.quality = fit.quality, adequacy = fit.adequacy, 
        predictivity = fit.predictivity, predictions = uit$predictions)
}
