PCAbipl.cat <-
function (X = Ocotea.data[, 3:8], G = NULL, X.new = NULL, scaled.mat = FALSE, 
    e.vects = 1:ncol(X), dim.biplot = c(2, 1, 3), adequacies.print = FALSE, 
    alpha = 0.95, alpha.3d = 0.7, aspect.3d = "iso", ax.col.3d = "black", 
    ax = 1:ncol(X), ax.name.col = rep("black", ncol(X)), ax.name.size = 0.65, 
    ax.type = c("predictive", "interpolative"), ax.col = list(ax.col = rep("grey", 
        ncol(X)), tickmarker.col = rep("grey", ncol(X)), marker.col = rep("black", 
        ncol(X))), between = c(1, -1, 0, 1), between.columns = -1, 
    cex.3d = 0.6, char.legend.size = c(1.2, 0.7), c.hull.n = 10, 
    colour.scheme = NULL, colours = c(1:8, 3:1), col.plane.3d = "lightgrey", 
    col.text.3d = "black", columns = 1, constant = 0.1, density.plot = FALSE, 
    exp.factor = 1.2, factor.x = 2, factor.y = 2, font.3d = 2, 
    ID.labs = FALSE, ID.3d = 1:nrow(X), label = TRUE, label.size = 0.6, 
    large.scale = FALSE, legend.type = c(means = FALSE, samples = FALSE, 
        bags = FALSE), line.length = c(1, 1), line.size = 2.5, 
    line.type = 1:ncol(G), line.width = 1, markers = TRUE, marker.size = 0.5, 
    max.num = 2500, means.plot = FALSE, n.int = rep(5, ncol(X)), 
    offset = rep(0, 4), ort.lty = 1, oblique.trans = NULL, orthog.transx = rep(0, 
        ncol(X)), orthog.transy = rep(0, ncol(X)), output = 1:9, 
    parlegendmar = c(3, 1, 3, 1), parplotmar = rep(3, 4), pch.means = 0:10, 
    pch.means.size = 1, pch.new = 1, pch.new.cols = "black", 
    pch.new.labels = NULL, pch.samples = 0:10, pch.samples.size = 1, 
    pos = c("Orthog", "Hor", "Paral"), predictions.3D = TRUE, 
    predictions.mean = NULL, predictions.sample = NULL, predictivity.print = FALSE, 
    quality.print = FALSE, reflect = c(FALSE, "x", "y"), rotate.degrees = 0, 
    select.origin = FALSE, size.ax.3d = 0.5, size.means.3d = 10, 
    size.points.3d = 5, specify.xaxis = NULL, specify.bags = NULL, 
    specify.classes = dimnames(G)[[2]], text.width.mult = 1, 
    Title = "", Titles.3d = c("", "", "x", "y", "z"), Tukey.median = TRUE) 
{
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    ax.type <- ax.type[1]
    pos <- pos[1]
    dim.biplot <- dim.biplot[1]
    if (is.null(G)) {
        G <- matrix(indmat(rep(1, nrow(X))), ncol = 1)
        dimnames(G) <- list(1:nrow(X), "AllData")
    }
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
    if (dim.biplot == 1) {
        out.list <- return(PCAbipl.OneDim(X = X, G = G, X.new = X.new, 
            e.vects = e.vects[1], adequacies.print = adequacies.print, 
            alpha = alpha, ax = ax, ax.col = ax.col, ax.name.col = ax.name.col, 
            ax.name.size = ax.name.size, ax.type = ax.type, between = between, 
            between.columns = between.columns, char.legend.size = char.legend.size, 
            c.hull.n = c.hull.n, colour.scheme = colour.scheme, 
            colours = colours, columns = columns, constant = constant, 
            density.plot = density.plot, exp.factor = exp.factor, 
            label = label, label.size = label.size, large.scale = large.scale, 
            legend.type = legend.type, line.length = line.length, 
            line.size = line.size, line.type = line.type, line.width = line.width, 
            markers = markers, marker.size = marker.size, max.num = max.num, 
            means.plot = means.plot, n.int = n.int, ort.lty = ort.lty, 
            parlegendmar = parlegendmar, parplotmar = parplotmar, 
            pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples = pch.samples, pch.samples.size = pch.samples.size, 
            pos = pos, predictions.mean = predictions.mean, predictions.sample = predictions.sample, 
            predictivity.print = predictivity.print, quality.print = quality.print, 
            offset = offset, reflect = reflect.mat, specify.bags = specify.bags, 
            specify.classes = specify.classes, Title = Title, 
            Tukey.median = Tukey.median))
    }
    if (dim.biplot == 3) {
        out.list <- return(PCAbipl.3d(X = X, G = G, X.new = X.new, 
            e.vects = e.vects[1:3], alpha.3d = alpha.3d, ax.col = ax.col[[1]], 
            aspect.3d = aspect.3d, ax.col.3d = ax.col.3d, cex.3d = cex.3d, 
            colour.scheme = colour.scheme, colours = c(4:12, 
                3:1), col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
            factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
            ID.labs = ID.labs, ID.3d = ID.3d, n.int = rep(5, 
                ncol(X)), plot.class.means = means.plot, predictions = predictions.3D, 
            size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
            size.means.3d = size.means.3d, specify.classes = specify.classes, 
            Titles.3d = Titles.3d))
    }
    if (dim.biplot == 2) {
        old.par <- par(no.readonly = TRUE)
        if (adequacies.print & predictivity.print) 
            stop("adequacies.print and predictivity.print cannot both be set to True")
        if (!is.null(G)) 
            if ((length(specify.bags) > ncol(G)) | (length(specify.classes) > 
                ncol(G))) 
                return(cat("Number of specified bags or specified classes must not be larger than the number of different classes\n"))
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
                dimnames(G) <- list(NULL, paste("class", 1:J, 
                  sep = ""))
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
        Vr.before.rotate <- svd.out$u[, e.vects[1:2]]
        eigval <- svd.out$d
        lambda.mat <- diag(eigval)
        eigval.r <- eigval[e.vects[1:2]]
        lambda.r.mat <- diag(eigval.r)
        fit.predictivity.mat <- diag(diag(Vr.before.rotate %*% 
            lambda.r.mat %*% t(Vr.before.rotate))) %*% solve(diag(diag(V.mat %*% 
            lambda.mat %*% t(V.mat))))
        fit.predictivity <- round(diag(fit.predictivity.mat), 
            digits = 3)
        names(fit.predictivity) <- dimnames(X)[[2]]
        fit.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
            eigval[e.vects[2]])/sum(eigval)) * 100, digits = 2), 
            "%")
        fit.adequacy <- round(diag(Vr.before.rotate %*% t(Vr.before.rotate)), 
            digits = 3)
        names(fit.adequacy) <- dimnames(X)[[2]]
        if (!is.null(specify.xaxis)) {
            if (!is.numeric(specify.xaxis)) 
                if (specify.xaxis == "maxpred") {
                  specify.xaxis <- (names(fit.predictivity))[fit.predictivity == 
                    max(fit.predictivity)]
                  specify.xaxis <- match(specify.xaxis, dimnames(X)[[2]])
                }
                else specify.xaxis <- match(specify.xaxis, dimnames(X)[[2]])
            radns <- -atan2(V.mat[specify.xaxis, e.vects[2]], 
                V.mat[specify.xaxis, e.vects[1]])
            rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
                cos(radns)), ncol = 2)
        }
        Vr <- Vr.before.rotate %*% rotate.mat %*% reflect.mat
        if (adequacies.print) 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", 
                fit.adequacy, ")", sep = "")
        if (predictivity.print) 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], " (", 
                round(fit.predictivity, digits = 2), ")", sep = "")
        if (!is.null(G)) {
            class.means.mat <- as.matrix(solve(t(G) %*% G) %*% 
                t(G) %*% unscaled.X, ncol = ncol(unscaled.X))
        }
        Z <- X %*% Vr
        if (!is.null(G)) {
            Z.means.mat <- scale(class.means.mat, means, sd) %*% 
                Vr
            Z.means.mat <- data.frame(Z.means.mat, pch.means = pch.means[1:J], 
                colr = as.character(colours[1:J]), line.type = line.type[1:J], 
                pch.means.size = pch.means.size, stringsAsFactors = FALSE)
            classnames <- dimnames(Z.means.mat)[[1]]
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
            if (ax.type == "interpolative") {
                if (any(orthog.transx != 0 | orthog.transy != 
                  0)) 
                  stop("Orthogonal translation only with prediction axes. \n")
                axes.rows <- Vr
            }
            else stop("ax.type is either 'interpolative' or 'predictive' (the default)")
        }
        z.axes <- lapply(1:p, function(j, unscaled.X, means, 
            sd, axes.rows, orthog.transx, orthog.transy) {
            phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                axes.rows %*% c(orthog.transx[j], orthog.transy[j])
            std.markers <- unique(sort(unscaled.X[, j]))
            interval <- (std.markers - means[j])/sd[j]
            axis.vals <- std.markers
            axis.points <- matrix(0, nrow = length(std.markers), 
                ncol = 4)
            axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 1]
            axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 2]
            if (!is.null(oblique.trans) & ax.type == "interpolative") {
                axis.points[, 1] <- (axis.vals) * axes.rows[j, 
                  1] - ((oblique.trans[j] - means[j])/sd[j]) * 
                  axes.rows[j, 1] + (((oblique.trans - means)/sd) %*% 
                  Vr[, 1])/p
                axis.points[, 2] <- (axis.vals) * axes.rows[j, 
                  2] - ((oblique.trans[j] - means[j])/sd[j]) * 
                  axes.rows[j, 2] + (((oblique.trans - means)/sd) %*% 
                  Vr[, 2])/p
            }
            axis.points[, 3] <- axis.vals * sd[j] + means[j]
            axis.points[, 4] <- 1
            return(axis.points)
        }, unscaled.X = unscaled.X, means = means, sd = sd, axes.rows = axes.rows, 
            orthog.transx = orthog.transx, orthog.transy = orthog.transy)
        uit <- drawbipl.catPCA(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
            p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
            alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
            pch.samples.size = pch.samples.size, specify.bags = select.numeric.bags, 
            label = label, markers = markers, Title = Title, 
            means.plot = means.plot, large.scale = large.scale, 
            specify.classes = select.numeric.classes, Tukey.median = Tukey.median, 
            Z.means.mat = Z.means.mat, offset = offset, pos = pos, 
            strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
            marker.size = marker.size, label.size = label.size, 
            exp.factor = exp.factor, line.width = line.width, 
            class.vec = class.vec, predictions.sample = predictions.sample, 
            predictions.mean = predictions.mean, ort.lty = ort.lty)
        if (!is.null(oblique.trans) & ax.type == "interpolative") 
            points(0, 0, pch = "+", cex = 2)
        if (!(is.null(X.new))) {
            if (is.null(pch.new)) {
                pch.new <- paste("N", 1:nrow(X.new), sep = "")
                Z.new <- scale(X.new, means, sd) %*% Vr
                text(Z.new, labels = pch.new, cex = pch.samples.size, 
                  col = "black")
            }
            else {
                Z.new <- scale(X.new, means, sd) %*% Vr
                points(Z.new, pch = pch.new, cex = pch.samples.size, 
                  col = pch.new.cols)
                text(Z.new, labels = pch.new.labels, pos = 1, 
                  cex = 0.75)
            }
        }
        if (!is.null(uit$legend.lab.bags)) 
            dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
        if (select.origin) {
            if (ax.type != "predictive") 
                stop("Select.origin only available for predictive biplot axes \n")
            cat("Move cursor where you want origin and press left mouse button \n")
            flush.console()
            origin.pos <- locator(1)
            dev.off()
            z.axes <- lapply(1:p, function(j, unscaled.X, means, 
                sd, axes.rows, n.int, orthog.transx, orthog.transy) {
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
                interval <- (std.markers - means[j])/sd[j]
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
                    1] - ((oblique.trans[j] - means[j])/sd[j]) * 
                    axes.rows[j, 1] + (((oblique.trans - means)/sd) %*% 
                    Vr[, 1])/p
                  axis.points[, 2] <- (axis.vals) * axes.rows[j, 
                    2] - ((oblique.trans[j] - means[j])/sd[j]) * 
                    axes.rows[j, 2] + (((oblique.trans - means)/sd) %*% 
                    Vr[, 2])/p
                }
                axis.points[, 3] <- axis.vals * sd[j] + means[j]
                axis.points[, 4] <- 0
                for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
                  3] - std.markers) == 0)) 
                  axis.points[i, 4] <- 1
                return(axis.points)
            }, unscaled.X = unscaled.X, means = means, sd = sd, 
                axes.rows = axes.rows, n.int = n.int, orthog.transx = rep(origin.pos$x, 
                  p), orthog.transy = rep(origin.pos$y, p))
            uit <- drawbipl.catPCA(Z, z.axes, z.axes.names = dimnames(X)[[2]], 
                p = p, ax = ax, ax.col = ax.col, ax.name.size = ax.name.size, 
                alpha = alpha, pch.means = pch.means, pch.means.size = pch.means.size, 
                pch.samples.size = pch.samples.size, specify.bags = select.numeric.bags, 
                label = label, markers = markers, Title = Title, 
                means.plot = means.plot, large.scale = large.scale, 
                specify.classes = select.numeric.classes, Tukey.median = Tukey.median, 
                Z.means.mat = Z.means.mat, offset = offset, pos = pos, 
                strepie = line.length, max.num = max.num, c.hull.n = c.hull.n, 
                marker.size = marker.size, label.size = label.size, 
                exp.factor = exp.factor, line.width = line.width, 
                class.vec = class.vec, predictions.sample = predictions.sample, 
                predictions.mean = predictions.mean, ort.lty = ort.lty)
            if (!is.null(oblique.trans) & ax.type == "interpolative") 
                points(0, 0, pch = "+", cex = 2)
            if (!(is.null(X.new))) {
                if (is.null(pch.new)) {
                  pch.new <- paste("N", 1:nrow(X.new), sep = "")
                  Z.new <- scale(X.new, means, sd) %*% Vr
                  text(Z.new, labels = pch.new, cex = pch.samples.size, 
                    col = "black")
                }
                else {
                  Z.new <- scale(X.new, means, sd) %*% Vr
                  points(Z.new, pch = pch.new, cex = pch.samples.size, 
                    col = pch.new.cols)
                }
            }
            if (!is.null(uit$legend.lab.bags)) 
                dimnames(G)[[2]][select.numeric.bags] <- uit$legend.lab.bags
        }
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
        out.list <- list(Z = Z, Z.axes = z.axes, V = V.mat, Eigenvectors = Vr, 
            e.vals = eigval, PCA.quality = fit.quality, adequacy = fit.adequacy, 
            predictivity = fit.predictivity, predictions = uit$predictions)
    }
    else stop("Argument dim.biplot must be specified to be 1, 2 or 3. \n")
    out.list[output]
}
