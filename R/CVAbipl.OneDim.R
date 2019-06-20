CVAbipl.OneDim <-
function (X = scale(Ocotea2.data[, 3:8]), G = indmat(Ocotea2.data[, 
    9]), X.new.samples = NULL, weightedCVA = c("weighted", "unweightI", 
    "unweightCent"), adequacies.print = FALSE, alpha = 0.5, ax = 1:ncol(X), 
    ax.name.col = rep("black", ncol(X)), ax.name.size = 0.65, 
    ax.col = list(ax.col = rep(2, ncol(X)), tickmarker.col = rep(3, 
        ncol(X)), marker.col = rep(3, ncol(X))), ax.type = c("predictive", 
        "interpolative"), between = c(1, -1, 0, 1), between.columns = -1, 
    c.hull.n = 10, columns = 1, char.legend.size = c(1.2, 0.7), 
    colours = c(4:12, 3:1), colours.means = c(4:12, 3:1), colour.scheme = NULL, 
    conf.alpha = NULL, constant = 0.1, density.plot = FALSE, 
    exp.factor = 1.2, e.vects = 1, label = TRUE, label.size = 0.6, 
    large.scale = FALSE, legend.type = c(means = TRUE, samples = FALSE, 
        bags = FALSE), line.size = 2.5, line.type = 1:ncol(G), 
    line.width = 1, line.length = c(1, 1), markers = TRUE, marker.size = 0.5, 
    max.num = 2500, means.plot = TRUE, n.int = rep(5, ncol(X)), 
    offset = rep(0.5, 4), offset.m = rep(0, ncol(X)), ort.lty = 1, 
    parlegendmar = c(3, 1, 3, 1), parplotmar = rep(3, 4), pch.means = 0:10, 
    pch.means.size = 1, pch.new = 1, pch.new.cols = "black", 
    pch.new.labels = NULL, pch.samples = 0:10, pch.samples.size = 1, 
    pos = c("Orthog", "Hor", "Paral"), pos.m = rep(1, ncol(X)), 
    side.label = rep("right", ncol(X)), specify.xaxis = FALSE, 
    predictions.mean = NULL, predictions.sample = NULL, predictivity.print = FALSE, 
    quality.print = FALSE, reflect = diag(1), specify.bags = NULL, 
    specify.classes = dimnames(G)[[2]], specify.ellipses = NULL, 
    text.width.mult = 1, Title = "", Tukey.median = FALSE) 
{
    dim.biplot <- 1
    density.plot <- density.plot[1]
    weightedCVA <- weightedCVA[1]
    X <- as.matrix(X)
    unscaled.X <- X
    X.cent <- scale(X, center = TRUE, scale = FALSE)
    X.means <- apply(X, 2, mean)
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    K <- min(p, J - 1)
    Gmat <- G
    Nmat <- t(Gmat) %*% Gmat
    XcentBar.groups <- solve(Nmat) %*% t(Gmat) %*% X.cent
    SSP.T <- t(X.cent) %*% X.cent
    SSP.B <- t(XcentBar.groups) %*% Nmat %*% XcentBar.groups
    SSP.W <- SSP.T - SSP.B
    Wmat <- SSP.W
    svd.Wmat <- svd(Wmat)
    lambdamatI <- diag(svd.Wmat$d)
    Lmat <- svd.Wmat$u %*% solve(sqrt(lambdamatI))
    if (weightedCVA == "weighted") {
        Cmat <- Nmat
        Csqrt <- sqrt(Nmat)
    }
    if (weightedCVA == "unweightedI") {
        Cmat <- diag(J)
        Csqrt <- Cmat
    }
    if (weightedCVA == "unweightedCent") {
        Cmat <- diag(J) - matrix(1/J, nrow = J, ncol = J)
        Csqrt <- Cmat
    }
    if (is.na(match(weightedCVA, c("weighted", "unweightedI", 
        "unweightedCent")))) 
        stop(" Argument 'weightedCVA' must be one of 'weighted','unweightedI','unweightedCent' ")
    svd.step2 <- svd(t(Lmat) %*% t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups %*% Lmat)
    Vmat <- svd.step2$v
    lambdamat <- diag(svd.step2$d)
    svd.2sided <- Eigen.twosided(t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups, Wmat)
    Mmat <- svd.2sided$W
    lambdamat.2sided <- svd.2sided$Lambda.mat
    vec.temp <- rep(0, p)
    vec.temp[1:dim.biplot] <- 1
    Jmat <- diag(vec.temp)
    XcentBarLVJ <- XcentBar.groups %*% Mmat %*% Jmat
    XcentLVJ <- X.cent %*% Mmat %*% Jmat
    XcentBarHat <- XcentBar.groups %*% Mmat %*% Jmat %*% solve(Mmat)
    AxisPredictivity <- diag(t(XcentBarHat) %*% Cmat %*% XcentBarHat)/diag(t(XcentBar.groups) %*% 
        Cmat %*% XcentBar.groups)
    ClassPredictivity <- diag(Csqrt %*% XcentBarHat %*% solve(Wmat) %*% 
        t(XcentBarHat) %*% Csqrt)/diag(Csqrt %*% XcentBar.groups %*% 
        solve(Wmat) %*% t(XcentBar.groups) %*% Csqrt)
    X.within <- (diag(n) - Gmat %*% (solve(Nmat)) %*% t(Gmat)) %*% 
        X.cent
    X.within.hat <- X.within %*% Mmat %*% Jmat %*% solve(Mmat)
    Within.group.axis.predictivity <- diag(t(X.within.hat) %*% 
        X.within.hat)/(diag(t(X.within) %*% X.within))
    Within.group.sample.predictivity <- diag((X.within.hat) %*% 
        solve(Wmat) %*% t(X.within.hat))/(diag((X.within) %*% 
        solve(Wmat) %*% t(X.within)))
    if (is.null(G)) 
        stop("G cannot be set to NULL. Groups must be specified for a CVA biplot \n")
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    ax.type <- ax.type[1]
    dim.biplot <- dim.biplot[1]
    pos <- pos[1]
    reflect.mat <- diag(dim.biplot)
    if ((reflect == "x") || (reflect == "y")) 
        if (dim.biplot == 1 && reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    rotate.mat <- diag(dim.biplot)
    if (dim.biplot == 1) {
        old.par <- par(no.readonly = TRUE)
        if (!is.null(predictions.sample) & !is.null(predictions.mean)) 
            stop("Either predictions.mean or predictions.sample or both must be NULL")
        if (adequacies.print & predictivity.print) 
            stop("adequacies.print and predictivity.print cannot both be set to TRUE")
        if (is.null(G)) {
            cat("Use function PCAbipl.boek for a PCA biplot or provide G for a CVA biplot\n")
            return(invisible())
        }
        if (!is.null(G)) 
            if ((length(specify.bags) > ncol(G)) | (length(specify.classes) > 
                ncol(G))) 
                stop("Number of specified bags or classes must not be larger than the number of different classes")
        if (ncol(G) == 2) 
            warning("Two class CVA biplot reduces to a line; second axis not uniquely defined, \n")
        par(pty = "m", mar = parplotmar)
        on.exit(par(old.par))
        if (is.null(dimnames(G))) 
            dimnames(G) <- list(NULL, paste("class", 1:J, sep = ""))
        if (length(dimnames(G)[[2]]) == 0) 
            dimnames(G)[[2]] <- paste("class", 1:J, sep = "")
        n.groups <- apply(G, 2, sum)
        if (!is.null(specify.bags)) 
            stop("Alpha bags constructed only for  2D biplots. \n")
        if (!is.null(specify.ellipses)) 
            stop("Kappa ellipses constructed only for  2D biplots. \n")
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
        if (nrow(X) != nrow(G)) {
            stop("Rows of X and G unequal!!")
        }
        class.vec <- apply(t(apply(G, 1, function(x) x == max(x))), 
            1, function(s, G) dimnames(G)[[2]][s], G = G)
        if (J > length(pch.means)) 
            stop(paste("Increase size of pch.means argument"))
        if (J > length(pch.samples)) 
            stop(paste("Increase size of pch.samples argument"))
        if (J > length(colours)) 
            stop(paste("Increase size of colours argument"))
        if (is.null(dimnames(X))) 
            dimnames(X) <- list(paste("s", 1:n, sep = ""), paste("V", 
                1:p, sep = ""))
        if (length(dimnames(X)[[1]]) == 0) 
            dimnames(X)[[1]] <- paste("s", 1:n, sep = "")
        if (length(dimnames(X)[[2]]) == 0) 
            dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
        means.mat <- XcentBar.groups
        dimnames(means.mat)[[1]] <- dimnames(G)[[2]]
        Br <- Mmat[, e.vects[1], drop = FALSE]
        Brr <- matrix(solve(Mmat)[e.vects[1], ], nrow = 1)
        if (!is.null(specify.xaxis)) 
            warning("xaxis rotation only for 2D biplots. \n")
        Z <- X.cent %*% Mmat[, e.vects[1], drop = FALSE] %*% 
            rotate.mat %*% reflect.mat
        CVA.mean <- means.mat %*% Mmat[, e.vects[1], drop = FALSE] %*% 
            rotate.mat %*% reflect.mat
        Z.means.mat <- data.frame(CVA.mean, 0, pch.means = pch.means[1:J], 
            colr = colours.means[1:J], line.type = line.type[1:J], 
            pch.means.size = pch.means.size, stringsAsFactors = FALSE)
        classnames <- dimnames(Z.means.mat)[[1]]
        Z <- data.frame(Z, 0, pch.samples = pch.samples[1], colr = colours[1], 
            line.type = line.type[1], stringsAsFactors = FALSE)
        if (J > 0) 
            for (j in 1:J) {
                Z[G[, j] == 1, 3] <- pch.samples[j]
                Z[G[, j] == 1, 4] <- as.character(colours[j])
                Z[G[, j] == 1, 5] <- line.type[j]
            }
        if (is.null(specify.bags)) {
            line.type <- NULL
            Z <- Z[, 1:4]
        }
        if (ax.type == "predictive") 
            axes.rows <- solve(diag(diag(t(Brr) %*% Brr))) %*% 
                t(Brr) %*% rotate.mat %*% reflect.mat
        else {
            if (ax.type == "interpolative") 
                axes.rows <- Br %*% rotate.mat %*% reflect.mat
            else stop("ax.type is either 'interpolative' or 'predictive' (the default). \n")
        }
        z.axes <- lapply(1:p, function(j, unscaled.X, means, 
            axes.rows, n.int) {
            number.points <- 100
            std.markers <- pretty(unscaled.X[, j], n = n.int[j])
            std.range <- c(min(std.markers), max(std.markers))
            std.markers.min <- std.markers - (std.range[2] - 
                std.range[1])
            std.markers.max <- std.markers + (std.range[2] - 
                std.range[1])
            std.markers <- c(std.markers, std.markers.min, std.markers.max)
            interval <- std.markers - means[j]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            number.points <- length(axis.vals)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            axis.points[, 1] <- axis.vals * axes.rows[j, 1]
            axis.points[, 2] <- 0
            axis.points[, 3] <- axis.vals + means[j]
            axis.points[, 4] <- 0
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, means = X.means, axes.rows = axes.rows, 
            n.int = n.int)
        lambda.vec <- zapsmall(diag(lambdamat))
        CVA.quality.canvar <- sum(lambda.vec[e.vects[1:dim.biplot]])/sum(lambda.vec)
        CVA.quality.origvar <- sum(diag(t(XcentBarHat) %*% Cmat %*% 
            XcentBarHat))/sum(diag(t(XcentBar.groups) %*% Cmat %*% 
            XcentBar.groups))
        adequacy.p <- diag(Br %*% t(Br))/diag(Mmat %*% t(Mmat))
        names(adequacy.p) <- dimnames(X)[[2]]
        if (adequacies.print) 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", 
                round(adequacy.p, digits = 3), ")", sep = "")
        uit <- drawbipl.onedim.bagalpha(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.col = ax.name.col, 
            ax.name.size = ax.name.size, alpha = alpha, conf.alpha = conf.alpha, 
            pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples.size = pch.samples.size, specify.bags = NULL, 
            label = label, markers = markers, Title = Title, 
            means.plot = means.plot, large.scale = large.scale, 
            specify.classes = select.numeric.classes, Tukey.median = Tukey.median, 
            Z.means.mat = Z.means.mat, offset = offset, offset.m = offset.m, 
            ort.lty = ort.lty, pos = pos, pos.m = pos.m, side.label = side.label, 
            strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, constant = constant, 
            density.plot = density.plot)
        if (!(is.null(X.new.samples))) {
            Z.new <- matrix((scale(X.new.samples, X.means, scale = FALSE) %*% 
                Br %*% rotate.mat %*% reflect.mat), ncol = 1)
            Z.new <- cbind(Z.new, 0)
            if (is.null(pch.new)) {
                pch.new <- paste("N", 1:nrow(X.new.samples), 
                  sep = "")
                points(Z.new, pch = pch.new, cex = pch.samples.size, 
                  col = pch.new.cols)
                text(Z.new, labels = pch.new, cex = pch.samples.size, 
                  col = "black")
            }
            else {
                points(Z.new, pch = pch.new, cex = pch.samples.size, 
                  col = pch.new.cols)
                text(Z.new, labels = pch.new.labels, pos = 1, 
                  cex = 0.75, col = pch.new.cols)
            }
        }
        if (any(legend.type)) {
            dev.new()
            blegend.colchar(quality.print = quality.print, QualityOfDisplay = display.quality, 
                classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J]), 
                pch.samples = as.vector(pch.samples[1:J]), colours = as.vector(colours[1:J]), 
                colours.means = as.vector(colours.means[1:J]), 
                line.type = as.vector(line.type[1:J]), line.width = line.width, 
                pch.samples.size = char.legend.size, legend.type = legend.type, 
                between = between, columns = columns, parlegendmar = parlegendmar, 
                line.size = line.size, between.columns = between.columns)
        }
    }
    list(Mmat = Mmat, MmatJmat = Br, CVA.quality.canvar = CVA.quality.canvar, 
        CVA.quality.origvar = CVA.quality.origvar, adequacy = adequacy.p, 
        predictions = uit$predictions, Z.means.mat = Z.means.mat)
}
