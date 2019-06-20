CVAbipl.density <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    X.new.samples = NULL, X.new.vars = NULL, e.vects = 1:ncol(X), 
    weightedCVA = c("weighted", "unweightI", "unweightCent"), 
    dim.biplot = c(2, 1, 3), adequacies.print = FALSE, alpha = 0.95, 
    alpha.3d = 0.7, aspect.3d = "iso", ax.col.3d = "black", ax = 1:sum(c(ncol(X), 
        ncol(X.new.vars))), ax.pred.above = NULL, ax.name.col = rep("black", 
        sum(c(ncol(X), ncol(X.new.vars)))), ax.name.size = 0.65, 
    ax.type = c("predictive", "interpolative"), ax.col = list(ax.col = rep("grey", 
        sum(c(ncol(X), ncol(X.new.vars)))), tickmarker.col = rep("grey", 
        sum(c(ncol(X), ncol(X.new.vars)))), marker.col = rep("black", 
        sum(c(ncol(X), ncol(X.new.vars))))), bandwidth = NULL, 
    basket.n = 36, basket.beta = 0.4, between = c(1, -1, 0, 1), 
    between.columns = -1, cex.3d = 0.6, char.legend.size = c(1.2, 
        0.7), c.hull.n = 10, colours.density = c("antiquewhite2", 
        "yellow", "red"), colour.scheme = NULL, colours = c(1:8, 
        3:1), colours.means = c(1:8, 3:1), col.plane.3d = "lightgrey", 
    col.text.3d = "black", columns = 1, conf.alpha = NULL, constant = 0.1, 
    cuts.density = 50, dens.ax.cex = 0.6, dens.ax.tcl = -0.2, 
    dens.mgp = c(0, -0.25, 0), draw.densitycontours = FALSE, 
    density.plot = c(NULL, "groups", "all"), exp.factor = 1.2, 
    ellipse.kappa = NULL, ellipse.alpha = NULL, factor.x = 2, 
    factor.y = 2, font.3d = 2, ID.labs = FALSE, ID.3d = 1:nrow(X), 
    label = TRUE, label.size = 0.6, large.scale = FALSE, layout.heights = c(100, 
        10), legend.type = c(means = FALSE, samples = FALSE, 
        bags = FALSE), line.length = c(1, 1), line.size = 2.5, 
    line.type = 1:ncol(G), line.width = 1, markers = TRUE, marker.size = 0.5, 
    max.num = 2500, means.plot = TRUE, n.int = rep(5, sum(c(ncol(X), 
        ncol(X.new.vars)))), offset = rep(0, 4), offset.m = rep(0.001, 
        sum(c(ncol(X), ncol(X.new.vars)))), ort.lty = 1, oblique.trans = NULL, 
    orthog.transx = rep(0, sum(c(ncol(X), ncol(X.new.vars)))), 
    orthog.transy = rep(0, sum(c(ncol(X), ncol(X.new.vars)))), 
    output = 1:10, parlegendmar = c(3, 1, 3, 1), parlegendmar.dens = c(2, 
        5, 0, 5), parplotmar = rep(3, 4), pch.means = 0:10, pch.means.size = 1, 
    pch.new = 1, pch.new.cols = "black", pch.new.labels = NULL, 
    pch.samples = 0:10, pch.samples.size = 1, pos = c("Orthog", 
        "Hor", "Paral"), pos.m = rep(1, sum(c(ncol(X), ncol(X.new.vars)))), 
    predictions.3D = TRUE, predictions.mean = NULL, predictions.sample = NULL, 
    predictivity.print = FALSE, quality.print = FALSE, reflect = c(FALSE, 
        "x", "y"), rotate.degrees = 0, select.origin = FALSE, 
    side.label = rep("right", sum(c(ncol(X), ncol(X.new.vars)))), 
    size.ax.3d = 0.5, size.means.3d = 10, size.points.3d = 5, 
    smooth.n = 100, specify.xaxis = NULL, specify.bags = NULL, 
    specify.classes = dimnames(G)[[2]], specify.density.class = "allsamples", 
    specify.ellipses = dimnames(G)[[2]], specify.beta.baskets = NULL, 
    text.width.mult = 1, Title = "", Titles.3d = c("", "", "x", 
        "y", "z"), Tukey.median = TRUE) 
{
    dim.biplot <- dim.biplot[1]
    if (dim.biplot != 2) 
        stop("CVAbipl.density only for 2D CVA biplot. \n\n Use CVAbipl for adding densities to 1D CVA biplot. \n")
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
    lambda.vec <- zapsmall(diag(lambdamat))
    CVA.quality.canvar <- sum(lambda.vec[e.vects[1:dim.biplot]])/sum(lambda.vec)
    CVA.quality.origvar <- sum(diag(t(XcentBarHat) %*% Cmat %*% 
        XcentBarHat))/sum(diag(t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups))
    Br <- Mmat[, e.vects[1:dim.biplot]]
    Brr <- solve(Mmat)[e.vects[1:dim.biplot], ]
    adequacy.p <- diag(Br %*% t(Br))/diag(Mmat %*% t(Mmat))
    names(adequacy.p) <- dimnames(X)[[2]]
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
    if (dim.biplot == 2) {
        reflect <- reflect[1]
        radns <- pi * rotate.degrees/180
        if (reflect == FALSE) 
            reflect.mat <- diag(2)
        else {
            if (reflect == "x") 
                reflect.mat <- diag(c(1, -1))
            else {
                if (reflect == "y") 
                  reflect.mat <- diag(c(-1, 1))
                else stop("Argument reflect can only set to be NULL, x or y. \n")
            }
        }
    }
    else reflect.mat <- diag(dim.biplot)
    if ((reflect == "x") || (reflect == "y")) 
        if (dim.biplot == 1 && reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    if (dim.biplot == 2) 
        rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
            cos(radns)), ncol = 2)
    else rotate.mat <- diag(dim.biplot)
    if (dim.biplot == 2) {
        old.par <- par(no.readonly = TRUE)
        layout(matrix(1:2, ncol = 1), heights = layout.heights)
        par(pty = "s", mar = parplotmar, cex = 1)
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
        par(pty = "s", mar = parplotmar)
        if (is.null(dimnames(G))) 
            dimnames(G) <- list(NULL, paste("class", 1:J, sep = ""))
        if (length(dimnames(G)[[2]]) == 0) 
            dimnames(G)[[2]] <- paste("class", 1:J, sep = "")
        n.groups <- apply(G, 2, sum)
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
        if (is.numeric(specify.ellipses)) 
            specify.ellipses <- dimnames(G)[[2]][specify.ellipses]
        if (!is.null(specify.ellipses)) {
            if (!is.numeric(specify.ellipses)) {
                test <- match(specify.ellipses, unique(as.character(dimnames(G)[[2]])), 
                  nomatch = 0)
                if (sum(test != 0) < length(specify.ellipses)) 
                  stop(" One or more of specified ellipses non-existent.  Check spelling")
            }
        }
        if (!is.null(specify.ellipses)) {
            select.numeric.ellipses <- rep(0, length(specify.ellipses))
            for (k in 1:length(specify.ellipses)) select.numeric.ellipses[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                specify.ellipses[k]]
        }
        else {
            if (is.null(specify.ellipses)) 
                select.numeric.ellipses <- NULL
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
        if (is.null(specify.density.class) | (length(specify.density.class) > 
            1)) 
            stop("Exactly one density class must be specified for the plot")
        if (specify.density.class == "allsamples") 
            select.numeric.density.class <- -10
        else {
            if (is.numeric(specify.density.class)) 
                specify.density.class <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                  specify.density.class]
            else {
                if (!is.numeric(specify.density.class)) {
                  test <- match(specify.density.class, unique(as.character(dimnames(G)[[2]])), 
                    nomatch = 0)
                  if (sum(test != 0) < length(specify.density.class)) 
                    stop("The specified density class is non-existent.  Check spelling")
                }
                select.numeric.density.class <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
                  specify.density.class]
            }
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
        dimnames(G)[[2]] <- paste(dimnames(G)[[2]], "; n = ", 
            n.groups, sep = "")
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
        if (!is.null(specify.xaxis)) {
            if (!is.numeric(specify.xaxis)) 
                specify.xaxis <- match(specify.xaxis, dimnames(X)[[2]])
            if (ax.type == "predictive") 
                radns <- -atan2((solve(diag(diag(t(Brr) %*% Brr))) %*% 
                  t(Brr))[specify.xaxis, e.vects[2]], (solve(diag(diag(t(Brr) %*% 
                  Brr))) %*% t(Brr))[specify.xaxis, e.vects[1]])
            if (ax.type == "interpolative") 
                radns <- -atan2(Br[specify.xaxis, e.vects[2]], 
                  Br[specify.xaxis, e.vects[1]])
            rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
                cos(radns)), ncol = 2)
        }
        Z <- X.cent %*% Mmat[, e.vects[1:dim.biplot]] %*% rotate.mat %*% 
            reflect.mat
        CVA.mean <- means.mat %*% Mmat[, e.vects[1:dim.biplot]] %*% 
            rotate.mat %*% reflect.mat
        Z.means.mat <- data.frame(CVA.mean, pch.samples.mean = pch.means[1:J], 
            colr = colours[1:J], line.type = line.type[1:J], 
            pch.means.size = pch.means.size, stringsAsFactors = FALSE)
        classnames <- dimnames(Z.means.mat)[[1]]
        Z <- data.frame(Z, pch.samples.samp = pch.samples[1], 
            colr = colours[1], line.type = line.type[1], stringsAsFactors = FALSE)
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
            axes.rows, n.int, orthog.transx, orthog.transy) {
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
            interval <- std.markers - means[j]
            axis.vals <- seq(from = min(interval), to = max(interval), 
                length = number.points)
            axis.vals <- sort(unique(c(axis.vals, interval)))
            number.points <- length(axis.vals)
            axis.points <- matrix(0, nrow = number.points, ncol = 4)
            axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 1]
            axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 2]
            if (!is.null(oblique.trans) & ax.type == "interpolative") {
                axis.points[, 1] <- (axis.vals) * axes.rows[j, 
                  1] - ((oblique.trans[j] - means[j])) * axes.rows[j, 
                  1]
                +(((oblique.trans - means)) %*% Br[, 1])/p
                axis.points[, 2] <- (axis.vals) * axes.rows[j, 
                  2] - ((oblique.trans[j] - means[j])) * axes.rows[j, 
                  2]
                +(((oblique.trans - means)) %*% Br[, 2])/p
            }
            axis.points[, 3] <- axis.vals + means[j]
            axis.points[, 4] <- 0
            for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                3] - std.markers) == 0)) 
                axis.points[i, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, means = X.means, axes.rows = axes.rows, 
            n.int = n.int, orthog.transx = orthog.transx, orthog.transy = orthog.transy)
        display.quality <- c(paste("CVA.quality.canvar =", round(CVA.quality.canvar * 
            100, digits = 2), "%"), paste("CVA.quality.origvar =", 
            round(CVA.quality.origvar * 100, digits = 2), "%"))
        if (adequacies.print) 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", 
                round(adequacy.p, digits = 3), ")", sep = "")
        uit <- drawbipl.bagalpha.density(Z = Z, z.axes = z.axes, 
            z.axes.names = dimnames(X)[[2]], p = p, ax = ax, 
            ax.col = ax.col, ax.name.col = ax.name.col, ax.name.size = ax.name.size, 
            alpha = alpha, basket.beta = basket.beta, basket.n = basket.n, 
            ellipse.kappa = ellipse.kappa, ellipse.alpha = ellipse.alpha, 
            pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples.size = pch.samples.size, specify.bags = select.numeric.bags, 
            specify.ellipses = select.numeric.ellipses, label = label, 
            markers = markers, Title = Title, means.plot = means.plot, 
            large.scale = large.scale, specify.classes = select.numeric.classes, 
            specify.beta.baskets = select.numeric.beta.baskets, 
            specify.density.class = select.numeric.density.class, 
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
        if (!is.null(oblique.trans) & ax.type == "interpolative") 
            points(0, 0, pch = "+", cex = 2)
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
            Z.new <- scale(X.new.samples, X.means, scale = FALSE) %*% 
                Br %*% rotate.mat %*% reflect.mat
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
        if (!is.null(uit$legend.lab.bags)) 
            dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
        if (select.origin) {
            z.axes <- lapply(1:p, function(j, unscaled.X, means, 
                axes.rows, n.int, orthog.transx, orthog.transy) {
                phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                  axes.rows %*% c(orthog.transx[j], orthog.transy[j])
                number.points <- 100
                std.markers <- pretty(unscaled.X[, j], n = n.int[j])
                std.range <- c(min(std.markers), max(std.markers))
                std.markers.min <- std.markers - (std.range[2] - 
                  std.range[1])
                std.markers.max <- std.markers + (std.range[2] - 
                  std.range[1])
                std.markers <- c(std.markers, std.markers.min, 
                  std.markers.max)
                interval <- std.markers - means[j]
                axis.vals <- seq(from = min(interval), to = max(interval), 
                  length = number.points)
                axis.vals <- sort(unique(c(axis.vals, interval)))
                number.points <- length(axis.vals)
                axis.points <- matrix(0, nrow = number.points, 
                  ncol = 4)
                axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                  phi.vec[j]) * axes.rows[j, 1]
                axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                  phi.vec[j]) * axes.rows[j, 2]
                if (!is.null(oblique.trans) & ax.type == "interpolative") {
                  axis.points[, 1] <- (axis.vals) * axes.rows[j, 
                    1] - ((oblique.trans[j] - means[j])) * axes.rows[j, 
                    1]
                  +(((oblique.trans - means)) %*% Br[, 1])/p
                  axis.points[, 2] <- (axis.vals) * axes.rows[j, 
                    2] - ((oblique.trans[j] - means[j])) * axes.rows[j, 
                    2]
                  +(((oblique.trans - means)) %*% Br[, 2])/p
                }
                axis.points[, 3] <- axis.vals + means[j]
                axis.points[, 4] <- 0
                for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                  3] - std.markers) == 0)) 
                  axis.points[i, 4] <- 1
                return(axis.points)
            }, unscaled.X = unscaled.X, means = X.means, axes.rows = axes.rows, 
                n.int = n.int, orthog.transx = rep(origin.pos$x, 
                  p), orthog.transy = rep(origin.pos$y, p))
            lambda.vec <- zapsmall(diag(lambdamat))
            display.quality <- c(paste("CVA.quality.canvar =", 
                round(CVA.quality.canvar * 100, digits = 2), 
                "%"), paste("CVA.quality.origvar =", round(CVA.quality.origvar * 
                100, digits = 2), "%"))
            if (adequacies.print) 
                dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", 
                  round(adequacy.p, digits = 3), ")", sep = "")
            par(pty = "s", mar = parplotmar)
            uit <- drawbipl.bagalpha.density(Z = Z, z.axes = z.axes, 
                z.axes.names = dimnames(X)[[2]], p = p, ax = ax, 
                ax.col = ax.col, ax.name.col = ax.name.col, ax.name.size = ax.name.size, 
                alpha = alpha, basket.beta = basket.beta, basket.n = basket.n, 
                ellipse.kappa = ellipse.kappa, ellipse.alpha = ellipse.alpha, 
                pch.means = pch.means, pch.means.size = pch.means.size, 
                pch.samples.size = pch.samples.size, specify.bags = select.numeric.bags, 
                specify.ellipses = select.numeric.ellipses, specify.beta.baskets = select.numeric.beta.baskets, 
                label = label, markers = markers, Title = Title, 
                means.plot = means.plot, large.scale = large.scale, 
                specify.classes = select.numeric.classes, specify.density.class = select.numeric.density.class, 
                Tukey.median = Tukey.median, Z.means.mat = Z.means.mat, 
                offset = offset, offset.m = offset.m, pos = pos, 
                pos.m = pos.m, strepie = line.length, max.num = max.num, 
                c.hull.n = c.hull.n, marker.size = marker.size, 
                label.size = label.size, exp.factor = exp.factor, 
                line.width = line.width, class.vec = class.vec, 
                side.label = side.label, smooth.n = smooth.n, 
                cuts.density = cuts.density, colours.density = colours.density, 
                draw.densitycontours = draw.densitycontours, 
                parlegendmar.dens = parlegendmar.dens, dens.ax.cex = dens.ax.cex, 
                dens.ax.tcl = dens.ax.tcl, dens.mgp = dens.mgp, 
                bandwidth = bandwidth)
            if (!(is.null(X.new.samples))) {
                Z.new <- scale(X.new.samples, X.means, scale = FALSE) %*% 
                  Br %*% rotate.mat %*% reflect.mat
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
            par(pty = "m", mar = parlegendmar.dens)
            draw.density.rect(levels.rect = uit$levels.rect, 
                col.use = uit$col.use, dens.ax.cex = dens.ax.cex, 
                dens.mgp = dens.mgp, dens.ax.tcl = dens.ax.tcl)
            if (!(is.null(X.new.samples))) {
                Z.new <- scale(X.new.samples, X.means, scale = FALSE) %*% 
                  Br %*% rotate.mat %*% reflect.mat
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
            if (!is.null(uit$legend.lab.bags)) 
                dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
        }
        if (any(legend.type)) {
            dev.new()
            blegend.colchar(quality.print = quality.print, QualityOfDisplay = display.quality, 
                classes = dimnames(G)[[2]], pch.means = as.vector(pch.means[1:J]), 
                pch.samples = as.vector(pch.samples[1:J]), colours = as.vector(colours[1:J]), 
                colours.means = colours.means, line.type = as.vector(line.type[1:J]), 
                line.width = line.width, pch.samples.size = char.legend.size, 
                legend.type = legend.type, between = between, 
                columns = columns, parlegendmar = parlegendmar, 
                line.size = line.size, between.columns = between.columns, 
                text.width.mult = text.width.mult)
        }
        out.list <- list(SSP.B = SSP.B, SSP.W = SSP.W, Lmat = Lmat, 
            Mmat = Mmat, lambda = diag(zapsmall(lambda.vec)), 
            CVA.quality.canvar = CVA.quality.canvar, CVA.quality.origvar = CVA.quality.origvar, 
            adequacy = adequacy.p, predictions = uit$predictions, 
            Z.means.mat = Z.means.mat)
    }
    else stop("Specify dim.biplot as 1 or 2 or 3. \n")
    out.list[output]
}
