cabipl.G <-
function (X, G = NULL, X.new = NULL, e.vects = 1:ncol(X), dim.biplot = c(2, 
    1, 3), alpha.3d = 0.7, aspect.3d = "iso", ax.col.3d = "black", 
    ax = 1:ncol(X), ax.name.col = rep("black", ncol(X)), ax.name.size = 0.65, 
    ax.col = list(ax.col = rep(8, ncol(X)), tickmarker.col = rep(8, 
        ncol(X)), marker.col = rep(1, ncol(X))), between = c(1, 
        -1, 0, 1), between.columns = -1, cex.3d = 0.6, char.legend.size = c(1.2, 
        0.7), colour.scheme = NULL, colours = c(4:12, 3:1), col.plane.3d = "lightgrey", 
    col.text.3d = "black", columns = 1, constant = 0.1, exp.factor = 1.2, 
    factor.x = 2, factor.y = 2, font.3d = 2, ID.labs = FALSE, 
    ID.3d = 1:nrow(X), label = TRUE, label.size = 0.6, legend.type = c(means = FALSE, 
        samples = TRUE, bags = FALSE), line.length = c(1, 1), 
    line.size = 2.5, line.type = 1:ncol(G), line.width = 1, markers = TRUE, 
    marker.size = 0.5, means.plot = TRUE, n.int = rep(5, ncol(X)), 
    offset = rep(0.5, 4), ort.lty = 1, parlegendmar = c(3, 1, 
        3, 1), parplotmar = rep(3, 4), pch.means = 0:10, pch.means.size = 1, 
    pch.samples = 0:10, pch.samples.size = 1, pos = c("Orthog", 
        "Hor", "Paral"), predictions.mat = FALSE, predictions.sample = NULL, 
    quality.print = FALSE, size.ax.3d = 0.5, size.means.3d = 10, 
    size.points.3d = 5, specify.classes = dimnames(G)[[2]], text.width.mult = 1, 
    Title = "", Titles.3d = c("", "", "x", "y", "z")) 
{
    pos <- pos[1]
    dim.biplot <- dim.biplot[1]
    if (!(dim.biplot[1] == "1" | dim.biplot[1] == "2" | dim.biplot[1] == 
        "3")) 
        stop("Specify dim.biplot as either 1 or 2 or 3 \n")
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    if (!is.null(X.new)) {
        X.new <- matrix(X.new, ncol = ncol(X))
        temp1 <- X.new/sum(X.new)
        temp.r <- apply(temp1, 1, sum)
        X.new.profiles <- sweep(temp1, 1, temp.r, "/")
    }
    else X.new.profiles <- NULL
    dim.biplot <- dim.biplot[1]
    if (is.null(G)) {
        G <- matrix(indmat(rep(1, nrow(X))), ncol = 1)
        dimnames(G) <- list(1:nrow(X), "AllData")
    }
    calibrated.axes <- function(p.col, X, means, sd, axes.rows, 
        n.int, dim.biplot) {
        out.list <- vector("list", p.col)
        for (j in 1:p.col) {
            number.points <- 50
            std.markers <- pretty(X[, j], n = n.int[j])
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
            if (dim.biplot == 2) {
                axis.points <- matrix(0, nrow = number.points, 
                  ncol = 4)
                axis.points[, 1] <- axis.vals * axes.rows[j, 
                  1]
                axis.points[, 2] <- axis.vals * axes.rows[j, 
                  2]
                axis.points[, 3] <- axis.vals * sd[j] + means[j]
                axis.points[, 4] <- 0
                for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                  3] - std.markers) == 0)) 
                  axis.points[i, 4] <- 1
                out.list[[j]] <- axis.points
            }
            if (dim.biplot == 1) {
                axis.points <- matrix(0, nrow = number.points, 
                  ncol = 4)
                axis.points[, 1] <- axis.vals * axes.rows[j, 
                  1]
                axis.points[, 2] <- 0
                axis.points[, 3] <- axis.vals * sd[j] + means[j]
                axis.points[, 4] <- 0
                for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                  3] - std.markers) == 0)) 
                  axis.points[i, 4] <- 1
                out.list[[j]] <- axis.points
            }
            if (dim.biplot == 3) {
                axis.points <- matrix(0, nrow = number.points, 
                  ncol = 5)
                axis.points[, 1] <- axis.vals * axes.rows[j, 
                  1]
                axis.points[, 2] <- axis.vals * axes.rows[j, 
                  2]
                axis.points[, 3] <- axis.vals * axes.rows[j, 
                  3]
                axis.points[, 4] <- axis.vals * sd[j] + means[j]
                axis.points[, 5] <- 0
                for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                  4] - std.markers) == 0)) 
                  axis.points[i, 5] <- 1
                out.list[[j]] <- axis.points
            }
        }
        out.list
    }
    old.par <- par(no.readonly = TRUE)
    if (!is.null(G)) 
        if (length(specify.classes) > ncol(G)) 
            return(cat("Number of specified classes must not be larger than the \n                   number of different classes \n"))
    par(pty = "s", mar = parplotmar)
    on.exit(par(old.par))
    n.row <- nrow(X)
    p.col <- ncol(X)
    J.classes <- ncol(G)
    n.groups <- apply(G, 2, sum)
    if (J.classes > length(pch.means)) 
        stop(paste("Increase size of pch.means argument"))
    if (J.classes > length(pch.samples)) 
        stop(paste("Increase size of pch.samples argument"))
    if (J.classes > length(colours)) 
        stop(paste("Increase size of colours argument"))
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste(1:n.row), paste("Col", 1:p.col, 
            sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste(1:n.row)
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("Col", 1:p.col, sep = "")
    if (!is.null(G)) {
        if (nrow(G) != n.row) 
            stop("number of rows of X and G differ")
        if (is.null(dimnames(G))) 
            dimnames(G) <- list(NULL, paste("class", 1:J.classes, 
                sep = ""))
        if (length(dimnames(G)[[2]]) == 0) 
            dimnames(G)[[2]] <- paste("class", 1:J.classes, sep = "")
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
        dimnames(G)[[2]] <- paste(dimnames(G)[[2]], "; n.rows = ", 
            n.groups, sep = "")
        if (length(dimnames(G)[[2]]) == 1) 
            class.vec <- rep(dimnames(G)[[2]], nrow(X))
        else {
            class.vec <- apply(t(apply(G, 1, function(x) x == 
                max(x))), 1, function(s, G) dimnames(G)[[2]][s], 
                G = G)
        }
    }
    Pmat <- as.matrix(X/sum(X))
    r.mass <- apply(Pmat, 1, sum)
    c.mass <- apply(Pmat, 2, sum)
    Dr <- diag(r.mass)
    Dc <- diag(c.mass)
    R.profile.mat <- sweep(Pmat, 1, r.mass, "/")
    Smat <- sqrt(solve(Dr)) %*% (Pmat - matrix(r.mass, ncol = 1) %*% 
        matrix(c.mass, nrow = 1)) %*% sqrt(solve(Dc))
    svd.Smat <- svd(Smat)
    Phi.mat <- sqrt(solve(Dr)) %*% svd.Smat$u
    GH.mat <- sqrt(solve(Dc)) %*% svd.Smat$v
    F.mat <- Phi.mat %*% diag(svd.Smat$d)
    G.mat <- GH.mat %*% diag(svd.Smat$d)
    Inertia.principal <- svd.Smat$d^2
    `%-Inertia` <- round(100 * Inertia.principal/(sum(Inertia.principal)), 
        digits = 4)
    GH.scaled <- sweep(GH.mat, 1, sqrt(c.mass), "*")
    scale.1 <- 1/diag(GH.scaled[, e.vects[1]] %*% t(GH.scaled[, 
        e.vects[1]]))
    scale.2 <- 1/diag(GH.scaled[, e.vects[1:2]] %*% t(GH.scaled[, 
        e.vects[1:2]]))
    scale.3 <- 1/diag(GH.scaled[, e.vects[1:3]] %*% t(GH.scaled[, 
        e.vects[1:3]]))
    if (ncol(G) == 1) {
        Z.means.mat <- matrix(c.mass, nrow = 1)
        dimnames(Z.means.mat) <- list(dimnames(G)[[2]], dimnames(X)[[2]])
    }
    else {
        if (ncol(G) > 1) {
            Z.means.mat <- matrix(0, nrow = J.classes, ncol = p.col)
            for (i in 1:J.classes) Z.means.mat[i, ] <- apply(Pmat[G[, 
                i], ], 2, sum)
            dimnames(Z.means.mat) <- list(dimnames(G)[[2]], dimnames(X)[[2]])
        }
        else stop("G must be a matrix with at least one column. \n")
    }
    Z.means.mat <- as.data.frame(Z.means.mat)
    if (dim.biplot == 3) {
        fit.quality <- paste("Quality of 3D display =", round(((Inertia.principal[e.vects[1]] + 
            Inertia.principal[e.vects[2]] + Inertia.principal[e.vects[3]])/sum(Inertia.principal)) * 
            100, digits = 4), "%")
        Zmat <- data.frame(F.mat[, e.vects[1:3]], colr = as.character(colours[1]), 
            stringsAsFactors = FALSE)
        dimnames(Zmat)[[1]] <- dimnames(X)[[1]]
        Z.means.mat.plot <- data.frame((Z.means.mat[, 1:p.col, 
            drop = F] %*% GH.mat[, 1:3]), pch.means = pch.means[1:J.classes], 
            colr = as.character(colours[1:J.classes]), stringsAsFactors = FALSE)
        z.axes <- calibrated.axes(p.col = p.col, X = R.profile.mat, 
            dim.biplot = dim.biplot, means = c.mass, sd = sqrt(c.mass), 
            axes.rows = diag(scale.3) %*% GH.scaled, n.int = n.int)
        draw.cabipl.3d(Z = Zmat, z.axes = z.axes, z.axes.names = dimnames(X)[[2]], 
            X.new.profiles = X.new.profiles, e.vects = e.vects[1:3], 
            alpha.3d = alpha.3d, ax.col = ax.col[[1]], ax.col.3d = ax.col.3d, 
            aspect.3d = aspect.3d, cex.3d = cex.3d, col.plane.3d = col.plane.3d, 
            col.text.3d = col.text.3d, class.vec = class.vec, 
            factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
            GH.mat = GH.mat, ID.labs = ID.labs, ID.3d = ID.3d, 
            p = p.col, size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
            size.means.3d = size.means.3d, specify.classes = select.numeric.classes, 
            Titles.3d = Titles.3d, Z.means.mat = Z.means.mat.plot)
        predictions <- NULL
        if (predictions.mat) 
            predictions.mat <- ca.predictions.mat(Pmat = Pmat, 
                e.vects = e.vects[1:3])
        else predictions.mat <- NULL
    }
    if (dim.biplot == 1) {
        fit.quality <- paste("Quality of 1D display =", round((Inertia.principal[e.vects[1]]/sum(Inertia.principal)) * 
            100, digits = 4), "%")
        Zmat <- data.frame(F.mat[, e.vects[1], drop = F], 0, 
            pch.samp = pch.samples[1], colr = as.character(colours[1]), 
            stringsAsFactors = FALSE)
        dimnames(Zmat)[[1]] <- dimnames(X)[[1]]
        z.axes <- calibrated.axes(p.col = p.col, X = R.profile.mat, 
            dim.biplot = dim.biplot, means = c.mass, sd = sqrt(c.mass), 
            axes.rows = diag(scale.1) %*% GH.scaled, n.int = n.int)
        out <- drawbipl.onedim.bagalpha(Z = Zmat, z.axes = z.axes, 
            z.axes.names = dimnames(X)[[2]], p = p.col, ax = ax, 
            ax.col = ax.col, label = label, pch.samples.size = pch.samples.size, 
            constant = constant, ax.name.size = ax.name.size, 
            Z.means.mat = Z.means.mat, specify.bags = NULL, specify.classes = select.numeric.classes, 
            offset = offset, pos = pos, strepie = line.length, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, predictions.sample = predictions.sample, 
            Title = Title)
        if (!(is.null(X.new))) {
            pch.new <- paste("N", 1:nrow(X.new.profiles), sep = "")
            Z.new <- cbind(X.new.profiles %*% GH.mat[, 1, drop = F], 
                0)
            text(Z.new, labels = pch.new, cex = pch.samples.size, 
                col = "black")
        }
        if (is.null(predictions.sample)) 
            predictions <- NULL
        else predictions <- out$predictions
        if (predictions.mat) 
            predictions.mat <- ca.predictions.mat(Pmat = Pmat, 
                e.vects = e.vects[1])
        else predictions.mat <- NULL
        if (any(legend.type)) {
            dev.new()
            blegend.colchar(quality.print = quality.print, QualityOfDisplay = fit.quality, 
                classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J.classes]), 
                pch.samples = as.vector(pch.samples[1:J.classes]), 
                colours = as.vector(colours[1:J.classes]), line.type = as.vector(line.type[1:J.classes]), 
                line.width = line.width, pch.samples.size = char.legend.size, 
                legend.type = legend.type, between = between, 
                columns = columns, parlegendmar = parlegendmar, 
                line.size = line.size, between.columns = between.columns)
        }
    }
    if (dim.biplot == 2) {
        fit.quality <- paste("Quality of 2D display =", round(((Inertia.principal[e.vects[1]] + 
            Inertia.principal[e.vects[2]])/sum(Inertia.principal)) * 
            100, digits = 4), "%")
        Zmat <- data.frame(F.mat[, e.vects[1:2]], pch.sampl = pch.samples[1], 
            colr = as.character(colours[1]), stringsAsFactors = FALSE)
        dimnames(Zmat)[[1]] <- dimnames(X)[[1]]
        z.axes <- calibrated.axes(p.col = p.col, X = R.profile.mat, 
            dim.biplot = dim.biplot, means = c.mass, sd = sqrt(c.mass), 
            axes.rows = diag(scale.2) %*% GH.scaled, n.int = n.int)
        out <- drawbipl.bagalpha(Z = Zmat, z.axes = z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p.col, ax = ax, ax.col = ax.col, label = label, 
            pch.samples.size = pch.samples.size, ax.name.size = ax.name.size, 
            Z.means.mat = Z.means.mat, specify.bags = NULL, specify.classes = select.numeric.classes, 
            offset = offset, pos = pos, strepie = line.length, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, predictions.sample = predictions.sample, 
            ort.lty = 1, Title = Title)
        if (is.null(predictions.sample)) 
            predictions <- NULL
        else predictions <- out$predictions
        if (!(is.null(X.new))) {
            pch.new <- paste("N", 1:nrow(X.new.profiles), sep = "")
            Z.new <- X.new.profiles %*% GH.mat[, 1:2]
            text(Z.new, labels = pch.new, cex = pch.samples.size, 
                col = "black")
        }
        if (predictions.mat) 
            predictions.mat <- ca.predictions.mat(Pmat = Pmat, 
                e.vects = e.vects[1:2])
        else predictions.mat <- NULL
        if (any(legend.type)) {
            dev.new()
            blegend.colchar(quality.print = quality.print, QualityOfDisplay = fit.quality, 
                classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J.classes]), 
                pch.samples = as.vector(pch.samples[1:J.classes]), 
                colours = as.vector(colours[1:J.classes]), line.type = as.vector(line.type[1:J.classes]), 
                line.width = line.width, pch.samples.size = char.legend.size, 
                legend.type = legend.type, between = between, 
                columns = columns, parlegendmar = parlegendmar, 
                line.size = line.size, between.columns = between.columns)
        }
    }
    list(R.profile.mat = R.profile.mat, c.mass = c.mass, r.mass = r.mass, 
        Phi.mat = Phi.mat, GH.mat = GH.mat, F.mat = F.mat, G.mat = G.mat, 
        GH.scaled = GH.scaled, `%-Inertia` = `%-Inertia`, QualityOfFit = fit.quality, 
        predictions = predictions, predictions.mat = predictions.mat)
}
