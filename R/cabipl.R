cabipl <- 
function (X = as.matrix(SACrime08.data), X.new = NULL, e.vects = 1:ncol(X), 
    ca.variant = c("PearsonResA", "PearsonResB", "IndepDev", 
        "ConRatioA", "ConRatioB", "Chisq2Rows", "Chisq2Cols", 
        "Corr", "RowProfA", "RowProfB", "PCA"), alpha.3d = 0.7, 
    aspect.3d = "iso", ax = 1:ncol(X), adjust.3d = c(0.5, 0.5), 
    ax.name.size = 0.75, ax.name.col = "black", axis.col = "red", 
    ax.col.3d = "black", cex.3d = 0.6, col.plane.3d = "lightgrey", 
    col.text.3d = "black", col.points.col = "red", col.points.size = 1, 
    col.points.label.size = 0.8, col.points.text = TRUE, ConRatioMinOne = FALSE, 
    constant = 0.05, dim.biplot = c(2, 1, 3), exp.factor = 1.2, 
    E.mat.null = FALSE, factor.x = 2, factor.y = 2, font.3d = 2, 
    ID.labs = FALSE, ID.3d = 1:nrow(X), lambda = FALSE, line.length = c(1, 
        1), logCRat = FALSE, markers = TRUE, marker.size = 0.5, 
    marker.col = "darkgrey", n.int = rep(3, ncol(X)), offset = rep(0.5, 
        4), offset.m = rep(1, ncol(X)), ort.lty = 1, parplotmar = rep(3, 
        4), pch.row.points = 16, pch.col.points = 15, PearsonRes.scaled.markers = FALSE, 
    plot.col.points = TRUE, pos = c("Orthog", "Hor", "Paral"), 
    pos.m = rep(1, ncol(X)), predictions.3d = TRUE, predictions.sample = NULL, 
    predictions.X.new = FALSE, propshift = 0, reflect = c(FALSE, 
        "x", "y"), rotate.degrees = 0, row.points.col = "green", 
    row.points.size = 1, row.points.label.size = 0.8, RowProf.scaled.markers = FALSE, 
    samples.plot = TRUE, scaled.mat = FALSE, show.origin = FALSE, 
    side.label = rep("right", ncol(X)), size.ax.3d = 1, size.points.3d = 10, 
    text.pos = c(1, 1), tick.marker.col = "darkgrey", Title = "", 
    Titles.3d = c("", "", "x", "y", "z"), output = 1:17, X.new.pch = 1, 
    X.new.col = "blue", X.new.pch.cex = 1.5, X.new.labels = paste("np", 
        1:(length(X.new)/ncol(X)), sep = "."), X.new.labels.cex = 0.6, 
    zoomval = NULL, predictions.allsamples.onaxis = NULL) 
{
	if (!is.null(X.new)) 
        if (!(((is.vector(X.new)) | (is.matrix(X.new))) & is.numeric(X.new))) 
            stop("X.new must be a numeric vector or matrix. \n")
    X.new.pars <- list(X.new.pch, X.new.col, X.new.pch.cex, X.new.labels, 
        X.new.labels.cex)
    points.new <- NULL
    testpredictions.new <- NULL
    col.plot.coords <- NULL
    calibrations.list <- NULL
    dim.biplot <- dim.biplot[1]
    if (dim.biplot != 1) {
        constant <- 0
        propshift <- 0
    }
    if (!(dim.biplot == 1 | dim.biplot == 2 | dim.biplot == 3)) 
        stop("Argument dim.biplot must be set to 1 or 2 or 3 \n")
    ca.variant <- ca.variant[1]
    if (is.na(match(ca.variant, c("PearsonResA", "PearsonResB", 
        "IndepDev", "ConRatioA", "ConRatioB", "Chisq2Rows", "Chisq2Cols", 
        "Corr", "RowProfA", "RowProfB", "PCA")))) 
        stop("ca variant must be one of: PearsonResA, PearsonResB, IndepDev, ConRatioA, ConRatioB, Chisq2Rows, Chisq2Cols, Corr, RowProf, PCA \n")
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    predictions <- NULL
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
    if ((reflect == "x") | (reflect == "y")) 
        if (dim.biplot == 1 && reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    if (dim.biplot == 2) 
        rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
            cos(radns)), ncol = 2)
    else rotate.mat <- diag(dim.biplot)
    n <- sum(X)
    p <- nrow(X)
    q <- ncol(X)
    if (length(axis.col) == 1) 
        axis.col <- rep(axis.col, q)
    if (length(ax.name.col) == 1) 
        ax.name.col <- rep(ax.name.col, q)
    if (length(marker.col) == 1) 
        marker.col <- rep(marker.col, q)
    if (length(tick.marker.col) == 1) 
        tick.marker.col <- rep(tick.marker.col, q)
    R.mat <- diag(apply(X, 1, sum))
    C.mat <- diag(apply(X, 2, sum))
    E.mat <- R.mat %*% matrix(1, nrow = p, ncol = q) %*% C.mat/n
    RHalf <- sqrt(R.mat)
    CHalf <- sqrt(C.mat)
    RMinHalf <- solve(RHalf)
    CMinHalf <- solve(CHalf)
    dev.mat <- X - E.mat
    Contingency.mat.prop <- solve(R.mat) %*% dev.mat %*% solve(C.mat)
    Contingency.mat <- sum(X) * Contingency.mat.prop + 1
    dimnames(Contingency.mat) <- dimnames(X)
    weighted.dev.mat <- (sqrt(solve(R.mat)) %*% dev.mat %*% sqrt(solve(C.mat)))
    if (E.mat.null) 
        weighted.dev.mat <- (sqrt(solve(R.mat)) %*% X %*% sqrt(solve(C.mat)))
    svd.weighted.dev.mat <- svd(weighted.dev.mat)
    sing.values <- svd.weighted.dev.mat$d
    Quality <- sum((sing.values^2)[e.vects][1:dim.biplot])/sum((sing.values^2)) * 
        100
    if (E.mat.null) 
        Quality <- sum((sing.values^2)[e.vects][1:dim.biplot])/sum((sing.values[-1]^2)) * 
            100
    dimnames(weighted.dev.mat) <- dimnames(X)
    if (!is.null(X.new)) {
        if (is.vector(X.new)) {
            R.new <- sum(X.new)
            weighted.dev.mat.new <- (1/sqrt(R.new)) * (X.new - 
                R.new * diag(C.mat)/n) %*% CMinHalf
        }
        else {
            R.new <- diag(apply(X.new, 1, sum))
            weighted.dev.mat.new <- solve(sqrt(R.new)) %*% (X.new - 
                R.new %*% matrix(1, nrow = nrow(X.new), ncol = ncol(X)) %*% 
                  (C.mat)/n) %*% CMinHalf
        }
    }
    else weighted.dev.mat.new <- NULL
    if (ca.variant == "PCA") {
        out <- PCAbipl(X = weighted.dev.mat, scaled.mat = scaled.mat, 
            ax = ax, colours = c(row.points.col, 4:12, 3:1), 
            pch.samples.size = row.points.size, pch.sample = pch.row.points, 
            Title = Title, specify.bags = NULL, label.size = row.points.label.size, 
            markers = FALSE, n.int = n.int, predictions.sample = predictions.sample, 
			predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
            offset = offset, reflect = reflect, rotate.degrees = rotate.degrees)
        Quality <- out$PCA.quality
    }
    dev.new()
    options(pty = "s")
    par(mar = parplotmar)
    if (ca.variant == "PearsonResA") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        USigHalf <- (svd.weighted.dev.mat$u %*% diag(sqrt(svd.weighted.dev.mat$d)))[, 
            vect.scaffolding] %*% rotate.mat %*% reflect.mat
        VSigHalf <- (svd.weighted.dev.mat$v %*% diag(sqrt(svd.weighted.dev.mat$d)))[, 
            vect.scaffolding] %*% rotate.mat %*% reflect.mat
        VSigHalf.Xnew <- (svd.weighted.dev.mat$v %*% diag(sqrt(svd.weighted.dev.mat$d)))[, 
            vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(VSigHalf * VSigHalf)/(q * sum(USigHalf * 
                USigHalf))
            lam <- sqrt(sqrt(lam.4))
            USigHalf <- USigHalf * lam
            VSigHalf <- VSigHalf/lam
        }
        if (dim.biplot == 1) {
            USigHalf <- cbind(USigHalf, 0)
            VSigHalf <- cbind(VSigHalf, 0)
        }
        lims <- c(min(USigHalf, VSigHalf), max(USigHalf, VSigHalf))
        plot(rbind(USigHalf[, 1:2], VSigHalf[, 1:2]), asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
            xaxs = "i", yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(USigHalf[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- USigHalf[, 1:2]
        samples.coords <- USigHalf
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- VSigHalf
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            centered.weighted.dev.mat <- scale(weighted.dev.mat, 
                scale = FALSE, center = TRUE)
            while (nrow(calibrations.mat) < 2) {
                if (PearsonRes.scaled.markers) {
                  markers.vals <- pretty(range(centered.weighted.dev.mat[, 
                    i] * sqrt(sum(X)) + attr(centered.weighted.dev.mat, 
                    "scaled:center")[i] * sqrt(sum(X))), n = n.int[i])
                  transformed.markers.vals <- markers.vals/sqrt(sum(X))
                }
                else {
                  markers.vals <- pretty(range(centered.weighted.dev.mat[, 
                    i] + attr(centered.weighted.dev.mat, "scaled:center")[i]), 
                    n = n.int[i])
                  transformed.markers.vals <- markers.vals
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(VSigHalf[i, 
                    1:2]^2)) * VSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(VSigHalf[i, 
                    1:2]^2)) * VSigHalf[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(VSigHalf[i, 
                    1:3]^2)) * VSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(VSigHalf[i, 
                    1:3]^2)) * VSigHalf[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(VSigHalf[i, 
                    1:3]^2)) * VSigHalf[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) 
                points.new <- cbind(lam * weighted.dev.mat.new %*% 
                  VSigHalf.Xnew %*% rotate.mat %*% reflect.mat * 
                  1/svd.weighted.dev.mat$d[vect.scaffolding], 
                  0)
            else points.new <- lam * weighted.dev.mat.new %*% 
                VSigHalf.Xnew %*% solve(diag(svd.weighted.dev.mat$d[vect.scaffolding])) %*% 
                rotate.mat %*% reflect.mat
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (PearsonRes.scaled.markers) 
                  predictions <- ca.predictivities(data = X)$X.hat[[3]] * 
                    sqrt(sum(X))
                else predictions <- ca.predictivities(data = X)$X.hat[[3]]
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            if (PearsonRes.scaled.markers) 
                testpredictions.new <- points.new %*% t(VSigHalf) * 
                  sqrt(sum(X))
            else testpredictions.new <- points.new %*% t(VSigHalf)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "PearsonResB") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        USig <- (svd.weighted.dev.mat$u %*% diag(svd.weighted.dev.mat$d))[, 
            vect.scaffolding] %*% rotate.mat %*% reflect.mat
        V <- svd.weighted.dev.mat$v[, vect.scaffolding] %*% rotate.mat %*% 
            reflect.mat
        V.Xnew <- svd.weighted.dev.mat$v[, vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(V * V)/(q * sum(USig * USig))
            lam <- sqrt(sqrt(lam.4))
            USig <- USig * lam
            V <- V/lam
        }
        if (dim.biplot == 1) {
            USig <- cbind(USig, 0)
            V <- cbind(V, 0)
        }
        lims <- c(min(USig, V), max(USig, V))
        plot(rbind(USig[, 1:2], V[, 1:2]), asp = 1, xlim = lims * 
            exp.factor, ylim = lims * exp.factor, type = "n", 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(USig[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- USig[, 1:2]
        samples.coords <- USig
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- V
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            centered.weighted.dev.mat <- scale(weighted.dev.mat, 
                scale = FALSE, center = TRUE)
            while (nrow(calibrations.mat) < 2) {
                if (PearsonRes.scaled.markers) {
                  markers.vals <- pretty(range(centered.weighted.dev.mat[, 
                    i] * sqrt(sum(X)) + attr(centered.weighted.dev.mat, 
                    "scaled:center")[i] * sqrt(sum(X))), n = n.int[i])
                  transformed.markers.vals <- markers.vals/sqrt(sum(X))
                }
                else {
                  markers.vals <- pretty(range(centered.weighted.dev.mat[, 
                    i] + attr(centered.weighted.dev.mat, "scaled:center")[i]), 
                    n = n.int[i])
                  transformed.markers.vals <- markers.vals
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) 
                points.new <- cbind(lam * weighted.dev.mat.new %*% 
                  V.Xnew %*% rotate.mat %*% reflect.mat, 0)
            else points.new <- lam * weighted.dev.mat.new %*% 
                V.Xnew %*% rotate.mat %*% reflect.mat
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (PearsonRes.scaled.markers) 
                  predictions <- ca.predictivities(data = X)$X.hat[[3]] * 
                    sqrt(sum(X))
                else predictions <- ca.predictivities(data = X)$X.hat[[3]]
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            if (PearsonRes.scaled.markers) 
                testpredictions.new <- points.new %*% t(V) * 
                  sqrt(sum(X))
            else testpredictions.new <- points.new %*% t(V)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "IndepDev") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RHalfUSigHalf <- RHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CHalfVSigHalf <- CHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CHalfVSigHalf.Xnew <- CHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(CHalfVSigHalf * CHalfVSigHalf)/(q * 
                sum(RHalfUSigHalf * RHalfUSigHalf))
            lam <- sqrt(sqrt(lam.4))
            RHalfUSigHalf <- RHalfUSigHalf * lam
            CHalfVSigHalf <- CHalfVSigHalf/lam
        }
        if (dim.biplot == 1) {
            RHalfUSigHalf <- cbind(RHalfUSigHalf, 0)
            CHalfVSigHalf <- cbind(CHalfVSigHalf, 0)
        }
        lims <- c(min(RHalfUSigHalf, CHalfVSigHalf), max(RHalfUSigHalf, 
            CHalfVSigHalf))
        plot(rbind(RHalfUSigHalf[, 1:2], CHalfVSigHalf[, 1:2]), 
            asp = 1, xlim = lims * exp.factor, ylim = lims * 
                exp.factor, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RHalfUSigHalf[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RHalfUSigHalf[, 1:2]
        samples.coords <- RHalfUSigHalf
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CHalfVSigHalf
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                markers.vals <- pretty(range(dev.mat[, i]), n = n.int[i])
                transformed.markers.vals <- markers.vals
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:2]^2)) * CHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:2]^2)) * CHalfVSigHalf[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) {
                points.new <- cbind(lam * sqrt(R.new) %*% weighted.dev.mat.new %*% 
                  CMinHalf %*% CHalfVSigHalf.Xnew %*% rotate.mat %*% 
                  reflect.mat * 1/svd.weighted.dev.mat$d[vect.scaffolding], 
                  0)
            }
            else {
                points.new <- lam * sqrt(R.new) %*% weighted.dev.mat.new %*% 
                  CMinHalf %*% CHalfVSigHalf.Xnew %*% solve(diag(svd.weighted.dev.mat$d[vect.scaffolding])) %*% 
                  rotate.mat %*% reflect.mat
            }
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                predictions <- RHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CHalf
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            testpredictions.new <- points.new %*% t(CHalfVSigHalf)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "ConRatioA") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSigHalf <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfVSigHalf <- CMinHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfVSigHalf.Xnew <- CMinHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(CMinHalfVSigHalf * CMinHalfVSigHalf)/(q * 
                sum(RMinHalfUSigHalf * RMinHalfUSigHalf))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSigHalf <- RMinHalfUSigHalf * lam
            CMinHalfVSigHalf <- CMinHalfVSigHalf/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSigHalf <- cbind(RMinHalfUSigHalf, 0)
            CMinHalfVSigHalf <- cbind(CMinHalfVSigHalf, 0)
        }
        lims <- c(min(RMinHalfUSigHalf, CMinHalfVSigHalf), max(RMinHalfUSigHalf, 
            CMinHalfVSigHalf))
        plot(rbind(RMinHalfUSigHalf, CMinHalfVSigHalf), asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
            xaxs = "i", yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RMinHalfUSigHalf[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RMinHalfUSigHalf[, 1:2]
        samples.coords <- RMinHalfUSigHalf
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CMinHalfVSigHalf
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                if (ConRatioMinOne) {
                  markers.vals <- pretty(range(Contingency.mat.prop[, 
                    i] * (sum(X))), n = n.int[i])
                  transformed.markers.vals <- markers.vals/(sum(X))
                }
                else {
                  markers.vals <- pretty(range(Contingency.mat.prop[, 
                    i] * (sum(X)) + 1), n = n.int[i])
                  transformed.markers.vals <- (markers.vals - 
                    1)/(sum(X))
                }
                if (logCRat) {
                  markers.vals <- pretty(range(2 * log(Contingency.mat.prop[, 
                    i] * (sum(X)) + 1)), n = n.int[i])
                  transformed.markers.vals <- (exp(markers.vals/2) - 
                    1)/(sum(X))
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfVSigHalf[i, 
                    1:2]^2)) * CMinHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfVSigHalf[i, 
                    1:2]^2)) * CMinHalfVSigHalf[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfVSigHalf[i, 
                    1:3]^2)) * CMinHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfVSigHalf[i, 
                    1:3]^2)) * CMinHalfVSigHalf[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CMinHalfVSigHalf[i, 
                    1:3]^2)) * CMinHalfVSigHalf[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) {
                points.new <- cbind(lam * solve(sqrt(R.new)) %*% 
                  weighted.dev.mat.new %*% CHalf %*% CMinHalfVSigHalf.Xnew %*% 
                  rotate.mat %*% reflect.mat * 1/svd.weighted.dev.mat$d[vect.scaffolding], 
                  0)
            }
            else {
                points.new <- lam * solve(sqrt(R.new)) %*% weighted.dev.mat.new %*% 
                  CHalf %*% CMinHalfVSigHalf.Xnew %*% solve(diag(svd.weighted.dev.mat$d[vect.scaffolding])) %*% 
                  rotate.mat %*% reflect.mat
            }
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (ConRatioMinOne) 
                  predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CMinHalf * sum(X)
                else predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CMinHalf * sum(X) + 1
                if (logCRat) 
                  predictions <- 2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CMinHalf * sum(X) + 1)
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if ((dim.biplot <= 2) && (logCRat) && (!is.null(predictions.sample))) {
            dimnames.out <- dimnames(out)
            if (dim.biplot == 1) 
                out <- t(2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[1]] %*% 
                  CMinHalf * sum(X) + 1)[predictions.sample, 
                  ax])
            if (dim.biplot == 2) 
                out <- t(2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[2]] %*% 
                  CMinHalf * sum(X) + 1)[predictions.sample, 
                  ax])
            dimnames(out) <- dimnames.out
        }
        if (!is.null(X.new)) {
            if (ConRatioMinOne) 
                testpredictions.new <- points.new %*% t(CMinHalfVSigHalf) * 
                  sum(X)
            else testpredictions.new <- points.new %*% t(CMinHalfVSigHalf) * 
                sum(X) + 1
            if (logCRat) 
                testpredictions.new <- 2 * log(points.new %*% 
                  t(CMinHalfVSigHalf) * sum(X) + 1)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "ConRatioB") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSig <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfV <- CMinHalf %*% svd.weighted.dev.mat$v[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfV.Xnew <- CMinHalf %*% (svd.weighted.dev.mat$v)[, 
            vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(CMinHalfV * CMinHalfV)/(q * sum(RMinHalfUSig * 
                RMinHalfUSig))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSig <- RMinHalfUSig * lam
            CMinHalfV <- CMinHalfV/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSig <- cbind(RMinHalfUSig, 0)
            CMinHalfV <- cbind(CMinHalfV, 0)
        }
        lims <- c(min(RMinHalfUSig, CMinHalfV), max(RMinHalfUSig, 
            CMinHalfV))
        plot(rbind(RMinHalfUSig, CMinHalfV), asp = 1, xlim = lims * 
            exp.factor, ylim = lims * exp.factor, type = "n", 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RMinHalfUSig[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RMinHalfUSig[, 1:2]
        samples.coords <- RMinHalfUSig
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CMinHalfV
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                if (ConRatioMinOne) {
                  markers.vals <- pretty(range(Contingency.mat.prop[, 
                    i] * (sum(X))), n = n.int[i])
                  transformed.markers.vals <- markers.vals/(sum(X))
                }
                else {
                  markers.vals <- pretty(range(Contingency.mat.prop[, 
                    i] * (sum(X)) + 1), n = n.int[i])
                  transformed.markers.vals <- (markers.vals - 
                    1)/(sum(X))
                }
                if (logCRat) {
                  markers.vals <- pretty(range(2 * log(Contingency.mat.prop[, 
                    i] * (sum(X)) + 1)), n = n.int[i])
                  transformed.markers.vals <- (exp(markers.vals/2) - 
                    1)/(sum(X))
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfV[i, 
                    1:2]^2)) * CMinHalfV[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfV[i, 
                    1:2]^2)) * CMinHalfV[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfV[i, 
                    1:3]^2)) * CMinHalfV[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfV[i, 
                    1:3]^2)) * CMinHalfV[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CMinHalfV[i, 
                    1:3]^2)) * CMinHalfV[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) {
                points.new <- cbind(lam * solve(sqrt(R.new)) %*% 
                  weighted.dev.mat.new %*% CHalf %*% CMinHalfV.Xnew %*% 
                  rotate.mat %*% reflect.mat, 0)
            }
            else {
                points.new <- lam * solve(sqrt(R.new)) %*% weighted.dev.mat.new %*% 
                  CHalf %*% CMinHalfV.Xnew %*% rotate.mat %*% 
                  reflect.mat
            }
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (ConRatioMinOne) 
                  predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CMinHalf * sum(X)
                else predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CMinHalf * sum(X) + 1
                if (logCRat) 
                  predictions <- 2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CMinHalf * sum(X) + 1)
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if ((dim.biplot <= 2) && (logCRat) && (!is.null(predictions.sample))) {
            dimnames.out <- dimnames(out)
            if (dim.biplot == 1) 
                out <- t(2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[1]] %*% 
                  CMinHalf * sum(X) + 1)[predictions.sample, 
                  ax])
            if (dim.biplot == 2) 
                out <- t(2 * log(RMinHalf %*% ca.predictivities(data = X)$X.hat[[2]] %*% 
                  CMinHalf * sum(X) + 1)[predictions.sample, 
                  ax])
            dimnames(out) <- dimnames.out
        }
        if (!is.null(X.new)) {
            if (ConRatioMinOne) 
                testpredictions.new <- points.new %*% t(CMinHalfV) * 
                  sum(X)
            else testpredictions.new <- points.new %*% t(CMinHalfV) * 
                sum(X) + 1
            if (logCRat) 
                testpredictions.new <- 2 * log(points.new %*% 
                  t(CMinHalfV) * sum(X) + 1)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "Chisq2Rows") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSig <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        V <- svd.weighted.dev.mat$v[, vect.scaffolding] %*% rotate.mat %*% 
            reflect.mat
        if (lambda) {
            lam.4 <- p * sum(V * V)/(q * sum(RMinHalfUSig * RMinHalfUSig))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSig <- RMinHalfUSig * lam
            V <- V/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSig <- cbind(RMinHalfUSig, 0)
            V <- cbind(V, 0)
        }
        lims <- c(min(RMinHalfUSig, V), max(RMinHalfUSig, V))
        plot(rbind(RMinHalfUSig, V), asp = 1, xlim = lims * exp.factor, 
            ylim = lims * exp.factor, type = "n", xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
            main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RMinHalfUSig[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RMinHalfUSig[, 1:2]
        samples.coords <- RMinHalfUSig
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- V
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                markers.vals <- pretty(range(RMinHalf %*% weighted.dev.mat[, 
                  i]), n = n.int[i])
                transformed.markers.vals <- markers.vals
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(V[i, 
                    1:3]^2)) * V[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) 
                points.new <- cbind(lam^2 * solve(sqrt(R.new)) %*% 
                  weighted.dev.mat.new %*% V[, 1], 0)
            else points.new <- lam^2 * solve(sqrt(R.new)) %*% 
                weighted.dev.mat.new %*% V
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
            if (!is.null(zoomval)) {
                zoomval <- zoom(zoomval)
                plot(rbind(RMinHalfUSig, V), asp = 1, xlim = zoomval[1:2], 
                  ylim = zoomval[3:4], type = "n", xlab = "", 
                  ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                  yaxs = "i", main = Title)
                usr <- par("usr")
                plotshift <- propshift * (usr[4] - usr[3])
                if (dim.biplot == 1) 
                  row.plot.coords <- cbind(RMinHalfUSig[, 1, 
                    drop = FALSE], plotshift)
                else row.plot.coords <- RMinHalfUSig[, 1:2]
                samples.coords <- RMinHalfUSig
                points(row.plot.coords, pch = pch.row.points, 
                  col = row.points.col, cex = row.points.size)
                text(x = row.plot.coords, labels = dimnames(X)[[1]], 
                  pos = text.pos[1], cex = row.points.label.size, 
                  col = row.points.col)
                col.plot.coords <- V
                if (plot.col.points) {
                  points(col.plot.coords, pch = pch.col.points, 
                    col = col.points.col, cex = col.points.size, 
                    )
                  if (col.points.text) 
                    text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                      pos = text.pos[2], cex = col.points.label.size, 
                      col = col.points.col)
                }
                if (predictions.X.new) {
                  out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                    row.plot.coords.mat = row.plot.coords, ax = ax, 
                    axes.names = ax.names, ax.name.size = ax.name.size, 
                    axis.col = axis.col, ax.name.col = ax.name.col, 
                    constant = constant, line.length = line.length, 
                    markers = markers, marker.size = marker.size, 
                    marker.col = marker.col, offset = offset, 
                    offset.m = offset.m, ort.lty = ort.lty, pos = pos, 
                    pos.m = pos.m, propshift = propshift, predictions.sample = predictions.sample, 
			        predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                    side.label = side.label, tick.marker.col = tick.marker.col, 
                    q = q, X.new.points = points.new))
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                else {
                  out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                    row.plot.coords.mat = row.plot.coords, ax = ax, 
                    axes.names = ax.names, ax.name.size = ax.name.size, 
                    axis.col = axis.col, ax.name.col = ax.name.col, 
                    constant = constant, line.length = line.length, 
                    markers = markers, marker.size = marker.size, 
                    marker.col = marker.col, offset = offset, 
                    offset.m = offset.m, ort.lty = ort.lty, pos = pos, 
                    pos.m = pos.m, propshift = propshift, predictions.sample = predictions.sample, 
			        predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                    side.label = side.label, tick.marker.col = tick.marker.col, 
                    q = q))
                  if (!is.null(X.new)) {
                    points(points.new, pch = X.new.pars[[1]], 
                      col = X.new.pars[[2]], cex = X.new.pars[[3]])
                    text(points.new, pos = 3, labels = X.new.pars[[4]], 
                      cex = X.new.pars[[5]])
                  }
                  if (show.origin) 
                    points(0, 0, pch = 3, cex = 1)
                }
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]]
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            testpredictions.new <- points.new %*% t(V)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "Chisq2Cols") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        U <- svd.weighted.dev.mat$u[, vect.scaffolding] %*% rotate.mat %*% 
            reflect.mat
        CMinHalfVSig <- (solve(CHalf) %*% svd.weighted.dev.mat$v %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfVSig.Xnew <- (solve(CHalf) %*% svd.weighted.dev.mat$v %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(CMinHalfVSig * CMinHalfVSig)/(q * 
                sum(U * U))
            lam <- sqrt(sqrt(lam.4))
            U <- U * lam
            CMinHalfVSig <- CMinHalfVSig/lam
        }
        if (dim.biplot == 1) {
            U <- cbind(U, 0)
            CMinHalfVSig <- cbind(CMinHalfVSig, 0)
        }
        lims <- c(min(U, CMinHalfVSig), max(U, CMinHalfVSig))
        plot(rbind(U, CMinHalfVSig), asp = 1, xlim = lims * exp.factor, 
            ylim = lims * exp.factor, type = "n", xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
            main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(U[, 1, drop = FALSE], plotshift)
        else row.plot.coords <- U[, 1:2]
        samples.coords <- U
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CMinHalfVSig
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size, )
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                markers.vals <- pretty(range(weighted.dev.mat %*% 
                  CMinHalf[, i]), n = n.int[i])
                transformed.markers.vals <- markers.vals
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfVSig[i, 
                    1:2]^2)) * CMinHalfVSig[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfVSig[i, 
                    1:2]^2)) * CMinHalfVSig[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CMinHalfVSig[i, 
                    1:3]^2)) * CMinHalfVSig[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CMinHalfVSig[i, 
                    1:3]^2)) * CMinHalfVSig[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CMinHalfVSig[i, 
                    1:3]^2)) * CMinHalfVSig[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) 
                points.new <- cbind(lam * weighted.dev.mat.new %*% 
                  CHalf %*% CMinHalfVSig.Xnew %*% rotate.mat %*% 
                  reflect.mat/(svd.weighted.dev.mat$d[vect.scaffolding]^(2)), 
                  0)
            else points.new <- lam * weighted.dev.mat.new %*% 
                CHalf %*% CMinHalfVSig.Xnew %*% diag((svd.weighted.dev.mat$d[vect.scaffolding]^(-2))) %*% 
                rotate.mat %*% reflect.mat
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
            if (!is.null(zoomval)) {
                zoomval <- zoom(zoomval)
                plot(rbind(U, CMinHalfVSig), asp = 1, xlim = zoomval[1:2], 
                  ylim = zoomval[3:4], type = "n", xlab = "", 
                  ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                  yaxs = "i", main = Title)
                usr <- par("usr")
                plotshift <- propshift * (usr[4] - usr[3])
                if (dim.biplot == 1) 
                  row.plot.coords <- cbind(U[, 1, drop = FALSE], 
                    plotshift)
                else row.plot.coords <- U[, 1:2]
                samples.coords <- U
                points(row.plot.coords, pch = pch.row.points, 
                  col = row.points.col, cex = row.points.size)
                text(x = row.plot.coords, labels = dimnames(X)[[1]], 
                  pos = text.pos[1], cex = row.points.label.size, 
                  col = row.points.col)
                col.plot.coords <- CMinHalfVSig
                if (plot.col.points) {
                  points(col.plot.coords, pch = pch.col.points, 
                    col = col.points.col, cex = col.points.size, 
                    )
                  if (col.points.text) 
                    text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                      pos = text.pos[2], cex = col.points.label.size, 
                      col = col.points.col)
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
                for (i in 1:q/2) lines(col.plot.coords[c(i, i + 
                  q/2), 1:2], lty = 1, col = "black")
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                predictions <- ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CMinHalf
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            testpredictions.new <- points.new %*% t(CMinHalfVSig)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "Corr") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSig <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CMinHalfVSig <- (solve(CHalf) %*% svd.weighted.dev.mat$v %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        if (lambda) {
            lam.4 <- p * sum(CMinHalfVSig * CMinHalfVSig)/(q * 
                sum(RMinHalfUSig * RMinHalfUSig))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSig <- RMinHalfUSig * lam
            CMinHalfVSig <- CMinHalfVSig/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSig <- cbind(RMinHalfUSig, 0)
            CMinHalfVSig <- cbind(CMinHalfVSig, 0)
        }
        lims <- c(min(RMinHalfUSig, CMinHalfVSig), max(RMinHalfUSig, 
            CMinHalfVSig))
        plot(rbind(RMinHalfUSig, CMinHalfVSig), asp = 1, xlim = lims * 
            1.2, ylim = lims * 1.2, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", main = Title)
        points(RMinHalfUSig[, 1:2], pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = RMinHalfUSig[, 1:2], labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        points(CMinHalfVSig[, 1:2], pch = pch.col.points, col = col.points.col, 
            cex = col.points.size)
        text(x = CMinHalfVSig[, 1:2], labels = dimnames(X)[[2]], 
            pos = text.pos[2], cex = col.points.label.size, col = col.points.col)
        if (show.origin) 
            points(0, 0, pch = 3, cex = 1)
        if (dim.biplot == 3) {
            graphics.off()
            require(rgl)
            graph.coord <- rbind(RMinHalfUSig[, 1:3], CMinHalfVSig[, 
                1:3])
            open3d()
            text3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
            points3d(graph.coord, size = 0)
            aspect3d(aspect.3d)
            if (samples.plot) 
                points3d(RMinHalfUSig[, 1:3], size = size.points.3d, 
                  color = row.points.col)
            points3d(CMinHalfVSig[, 1:3], size = size.points.3d, 
                color = col.points.col)
            if (ID.labs) {
                text3d(RMinHalfUSig[, 1:3], text = ID.3d, colour = row.points.col, 
                  font = font.3d, cex = cex.3d)
                text3d(CMinHalfVSig[, 1:3], text = dimnames(X)[[2]], 
                  col = col.points.col, font = font.3d, cex = cex.3d)
            }
            rgl.surface(x = c(min(graph.coord[, 1]), max(graph.coord[, 
                1])) * factor.x, z = c(min(graph.coord[, 2]), 
                max(graph.coord[, 2])) * factor.y, y = matrix(0, 
                nrow = 2, ncol = 2), coords = c(1, 3, 2), col = col.plane.3d, 
                alpha = alpha.3d)
            axes3d(col = ax.col.3d)
            title3d(main = Titles.3d[1], sub = Titles.3d[2], 
                xlab = Titles.3d[3], ylab = Titles.3d[4], zlab = Titles.3d[5])
        }
        out <- NULL
    }
    if (ca.variant == "RowProfA") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSigHalf <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CHalfVSigHalf <- CHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CHalfVSigHalf.Xnew <- CHalf %*% (svd.weighted.dev.mat$v %*% 
            diag(sqrt(svd.weighted.dev.mat$d)))[, vect.scaffolding]
        if (lambda) {
            lam.4 <- p * sum(CHalfVSigHalf * CHalfVSigHalf)/(q * 
                sum(RMinHalfUSigHalf * RMinHalfUSigHalf))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSigHalf <- RMinHalfUSigHalf * lam
            CHalfVSigHalf <- CHalfVSigHalf/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSigHalf <- cbind(RMinHalfUSigHalf, 0)
            CHalfVSigHalf <- cbind(CHalfVSigHalf, 0)
        }
        lims <- c(min(RMinHalfUSigHalf, CHalfVSigHalf), max(RMinHalfUSigHalf, 
            CHalfVSigHalf))
        plot(rbind(RMinHalfUSigHalf, CHalfVSigHalf), asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
            xaxs = "i", yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RMinHalfUSigHalf[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RMinHalfUSigHalf[, 1:2]
        samples.coords <- RMinHalfUSigHalf
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CHalfVSigHalf
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size)
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                if (RowProf.scaled.markers) {
                  marginal <- (diag(C.mat)[i])/sum(X)
                  markers.vals <- pretty(range((solve(R.mat)) %*% 
                    dev.mat[, i] + marginal), n = n.int[i])
                  transformed.markers.vals <- markers.vals - 
                    marginal
                }
                else {
                  markers.vals <- pretty(range((solve(R.mat)) %*% 
                    dev.mat[, i]), n = n.int[i])
                  transformed.markers.vals <- markers.vals
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:2]^2)) * CHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:2]^2)) * CHalfVSigHalf[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CHalfVSigHalf[i, 
                    1:3]^2)) * CHalfVSigHalf[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            if (dim.biplot == 1) {
                points.new <- cbind(lam * solve(sqrt(R.new)) %*% 
                  weighted.dev.mat.new %*% solve(CHalf) %*% CHalfVSigHalf.Xnew %*% 
                  rotate.mat %*% reflect.mat/svd.weighted.dev.mat$d[vect.scaffolding], 
                  0)
            }
            else points.new <- {
                lam * solve(sqrt(R.new)) %*% weighted.dev.mat.new %*% 
                  solve(CHalf) %*% CHalfVSigHalf.Xnew %*% diag((svd.weighted.dev.mat$d[vect.scaffolding]^(-1))) %*% 
                  rotate.mat %*% reflect.mat
            }
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (RowProf.scaled.markers) 
                  predictions <- scale(RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CHalf, scale = FALSE, center = -diag(C.mat)/sum(X))
                else predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CHalf
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            if (RowProf.scaled.markers) 
                testpredictions.new <- scale(points.new %*% t(CHalfVSigHalf), 
                  scale = FALSE, center = -diag(C.mat)/sum(X))
            else testpredictions.new <- points.new %*% t(CHalfVSigHalf)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (ca.variant == "RowProfB") {
        vect.scaffolding <- e.vects[1:dim.biplot]
        Biplot.quality <- sum((sing.values^2)[vect.scaffolding])/sum((sing.values^2)) * 
            100
        RMinHalfUSig <- RMinHalf %*% (svd.weighted.dev.mat$u %*% 
            diag(svd.weighted.dev.mat$d))[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        CHalfV <- CHalf %*% svd.weighted.dev.mat$v[, vect.scaffolding] %*% 
            rotate.mat %*% reflect.mat
        if (lambda) {
            lam.4 <- p * sum(CHalfV * CHalfV)/(q * sum(RMinHalfUSig * 
                RMinHalfUSig))
            lam <- sqrt(sqrt(lam.4))
            RMinHalfUSig <- RMinHalfUSig * lam
            CHalfV <- CHalfV/lam
        }
        if (dim.biplot == 1) {
            RMinHalfUSig <- cbind(RMinHalfUSig, 0)
            CHalfV <- cbind(CHalfV, 0)
        }
        lims <- c(min(RMinHalfUSig, CHalfV), max(RMinHalfUSig, 
            CHalfV))
        plot(rbind(RMinHalfUSig, CHalfV), asp = 1, xlim = lims * 
            exp.factor, ylim = lims * exp.factor, type = "n", 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) 
            row.plot.coords <- cbind(RMinHalfUSig[, 1, drop = FALSE], 
                plotshift)
        else row.plot.coords <- RMinHalfUSig[, 1:2]
        samples.coords <- RMinHalfUSig
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        col.plot.coords <- CHalfV
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = col.points.col, 
                cex = col.points.size)
            if (col.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = col.points.label.size, 
                  col = col.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            while (nrow(calibrations.mat) < 2) {
                if (RowProf.scaled.markers) {
                  marginal <- (diag(C.mat)[i])/sum(X)
                  markers.vals <- pretty(range((solve(R.mat)) %*% 
                    dev.mat[, i] + marginal), n = n.int[i])
                  transformed.markers.vals <- markers.vals - 
                    marginal
                }
                else {
                  markers.vals <- pretty(range((solve(R.mat)) %*% 
                    dev.mat[, i]), n = n.int[i])
                  transformed.markers.vals <- markers.vals
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfV[i, 
                    1:2]^2)) * CHalfV[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfV[i, 
                    1:2]^2)) * CHalfV[i, 2]
                  calibrations.mat <- matrix(c(calibrations.x, 
                    calibrations.y, markers.vals), ncol = 3)
                  criterion.x <- calibrations.x > usr[1] & calibrations.x < 
                    usr[2]
                  criterion.y <- calibrations.y > usr[3] & calibrations.y < 
                    usr[4]
                  calibrations.mat <- calibrations.mat[criterion.x & 
                    criterion.y, , drop = FALSE]
                }
                else {
                  calibrations.x <- (transformed.markers.vals/sum(CHalfV[i, 
                    1:3]^2)) * CHalfV[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(CHalfV[i, 
                    1:3]^2)) * CHalfV[i, 2]
                  calibrations.z <- (transformed.markers.vals/sum(CHalfV[i, 
                    1:3]^2)) * CHalfV[i, 3]
                  calibrations.mat <- cbind(calibrations.x, calibrations.y, 
                    calibrations.z, markers.vals)
                }
                n.int[i] <- 2 * n.int[i]
            }
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (!is.null(X.new)) {
            if (lambda) 
                lam <- lam
            else lam <- 1
            points.new <- lam^2 * solve(sqrt(R.new)) %*% weighted.dev.mat.new %*% 
                solve(CHalf) %*% CHalfV
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            if (predictions.X.new) {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q, X.new.points = points.new))
                points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                  cex = X.new.pars[[3]])
                text(points.new, pos = 3, labels = X.new.pars[[4]], 
                  cex = X.new.pars[[5]])
            }
            else {
                out <- (drawbipl.ca(calibrations.list = calibrations.list, 
                  row.plot.coords.mat = row.plot.coords, ax = ax, 
                  axes.names = ax.names, ax.name.size = ax.name.size, 
                  axis.col = axis.col, ax.name.col = ax.name.col, 
                  constant = constant, line.length = line.length, 
                  markers = markers, marker.size = marker.size, 
                  marker.col = marker.col, offset = offset, offset.m = offset.m, 
                  ort.lty = ort.lty, pos = pos, pos.m = pos.m, 
                  propshift = propshift, predictions.sample = predictions.sample, 
			      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                  side.label = side.label, tick.marker.col = tick.marker.col, 
                  q = q))
                if (!is.null(X.new)) {
                  points(points.new, pch = X.new.pars[[1]], col = X.new.pars[[2]], 
                    cex = X.new.pars[[3]])
                  text(points.new, pos = 3, labels = X.new.pars[[4]], 
                    cex = X.new.pars[[5]])
                }
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
            }
        }
        if (dim.biplot == 3) {
            out <- (drawbipl.3dim.ca(Z = samples.coords, z.axes = calibrations.list, 
                X.new = points.new, X.new.pars = X.new.pars, 
                axes.names = ax.names, adjust.3 = adjust.3d, 
                alpha.3d = alpha.3d, ax.col = axis.col, ax.col.3d = ax.col.3d, 
                aspect.3d = aspect.3d, cex.3d = cex.3d, col.samples = row.points.col, 
                col.plane.3d = col.plane.3d, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3d) {
                if (RowProf.scaled.markers) 
                  predictions <- scale(RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                    CHalf, scale = FALSE, center = -diag(C.mat)/sum(X))
                else predictions <- RMinHalf %*% ca.predictivities(data = X)$X.hat[[3]] %*% 
                  CHalf
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
        if (!is.null(X.new)) {
            if (RowProf.scaled.markers) 
                testpredictions.new <- scale(points.new %*% t(CHalfV), 
                  scale = FALSE, center = -diag(C.mat)/sum(X))
            else testpredictions.new <- points.new %*% t(CHalfV)
            dimnames(testpredictions.new) <- list(X.new.pars[[4]], 
                dimnames(X)[[2]])
        }
    }
    if (!lambda) 
        lam <- NULL
    out.list <- list(out = out, X.mat = X, E.mat = E.mat, R.mat = R.mat, 
        C.mat = C.mat, dev.mat = dev.mat, weighted.dev.mat = weighted.dev.mat, 
        Contingency.mat.prop = Contingency.mat.prop, Contingency.mat = Contingency.mat, 
        svd.weighted.dev.mat = svd.weighted.dev.mat, lambda = lam, 
        Quality = Quality, predictions = predictions, weighted.dev.mat.new = weighted.dev.mat.new, 
        predictions.new = testpredictions.new, column.plot.coords = col.plot.coords, 
        calibrations.list = calibrations.list)
    out.list[output]
}
