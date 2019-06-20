biadbipl <-
function (X = wheat.data, X.new.rows = NULL, X.new.columns = NULL, 
    e.vects = 1:ncol(X), add.maineffect = FALSE, show.maineffect = FALSE, 
    show.origin = FALSE, biad.variant = c("InteractionMat", "XminMeanMat", 
        "XMat"), SigmaHalf = FALSE, circle.proj = NULL, expand.markers = NULL, 
    predictivity.print = FALSE, adj.3d = 0.5, alpha.3d = 0.7, 
    aspect.3d = "iso", ax = 1:ncol(X), ax.name.size = 0.75, ax.name.col = "black", 
    axis.col = "red", ax.col.3d = "black", cex.3d = 0.6, col.plane.3d = "lightgrey", 
    col.text.3d = "black", column.points.col = rep("red", ncol(X)), 
    column.points.size = 1, column.points.label.size = 0.8, column.points.text = TRUE, 
    constant = 0.05, dim.biplot = c(2, 1, 3), exp.factor = 1.2, 
    factor.x = 2, factor.y = 2, font.3d = 2, ID.labs = FALSE, 
    ID.3d = 1:nrow(X), lambda = FALSE, legend.show = FALSE, legend.columns = c(1, 
        1), line.length = c(1, 1), markers = TRUE, marker.size = 0.5, 
    marker.col = "grey", n.int = rep(3, ncol(X)), offset = rep(0.5, 
        4), offset.m = rep(0.5, sum(ncol(X), ncol(X.new.columns))), 
    ort.lty = 1, parplotmar = rep(3, 4), pch.row.points = 16, 
    pch.col.points = 15, plot.col.points = TRUE, pos = c("Orthog", 
        "Hor", "Paral"), pos.m = rep(1, sum(ncol(X), ncol(X.new.columns))), 
    predictions.3D = TRUE, predictions.sample = NULL, predictions.allsamples.onaxis = NULL, 
    propshift = 0, reflect = c(FALSE, "x", "y"), rotate.degrees = 0, 
    row.points.col = rep("green", nrow(X)), row.points.size = 1, 
    row.points.label.size = 0.8, select.origin = FALSE, samples.plot = TRUE, 
    side.label = rep("right", sum(ncol(X), ncol(X.new.columns))), 
    size.ax.3d = 1, size.points.3d = 10, text.pos = c(1, 1), 
    tick.marker.col = "grey", Title = "", Titles.3d = c("", "", 
        "x", "y", "z"), output = 1:12, X.new.columns.pch = 1, 
    X.new.columns.col = "blue", X.new.columns.pch.cex = 1.5, 
    X.new.columns.labels = dimnames(X.new.columns)[[2]], X.new.columns.labels.cex = 0.6, 
    X.new.rows.pch = 1, X.new.rows.col = "blue", X.new.rows.pch.cex = 1.5, 
    X.new.rows.labels = dimnames(X.new.rows)[[1]], X.new.rows.labels.cex = 0.6, 
    zoomval = NULL) 
{
    n.int.original <- n.int
    if (!is.null(X.new.rows)) 
        if (!(((is.vector(X.new.rows)) | (is.matrix(X.new.rows))) & 
            is.numeric(X.new.rows))) 
            stop("X.new.rows must be a numeric vector or matrix. \n")
    if (is.vector(X.new.rows) & (length(X.new.rows) != ncol(X))) 
        stop("X.new.rows not of correct size")
    if (is.vector(X.new.rows)) 
        X.new.rows <- matrix(X.new.rows, nrow = 1, ncol = ncol(X))
    X.new.rows.pars <- list(X.new.rows.pch, X.new.rows.col, X.new.rows.pch.cex, 
        X.new.rows.labels, X.new.rows.labels.cex)
    if (!is.null(X.new.columns)) 
        if (!(((is.vector(X.new.columns)) | (is.matrix(X.new.columns))) & 
            is.numeric(X.new.columns))) 
            stop("X.new.columns must be a numeric vector or matrix. \n")
    if (is.vector(X.new.columns) & (length(X.new.columns) != 
        nrow(X))) 
        stop("X.new.columns not of correct size")
    if (is.vector(X.new.columns)) 
        X.new.columns <- matrix(X.new.columns, ncol = 1, nrow = nrow(X))
    X.new.columns.pars <- list(X.new.columns.pch, X.new.columns.col, 
        X.new.columns.pch.cex, X.new.columns.labels, X.new.columns.labels.cex)
    points.new <- NULL
    testpredictions.new <- NULL
    dim.biplot <- dim.biplot[1]
    if (dim.biplot == 1) 
        circle.proj <- NULL
    if (dim.biplot != 1) {
        constant <- 0
        propshift <- 0
    }
    if (!(dim.biplot == 1 | dim.biplot == 2 | dim.biplot == 3)) 
        stop("Argument dim.biplot must be set to 1 or 2 or 3 \n")
    biad.variant <- biad.variant[1]
    if (is.na(match(biad.variant, c("Xmat", "XminMeanMat", "InteractionMat")))) 
        stop("biad.variant  must be one of: Xmat, XminMeanMat, InteractionMat\n")
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    predictions <- NULL
    reflect <- reflect[1]
    if (dim.biplot == 2) {
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
        if (dim.biplot == 1 & reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    if (dim.biplot == 2) 
        rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
            cos(radns)), ncol = 2)
    else rotate.mat <- diag(dim.biplot)
    Xmean <- mean(X)
    p <- nrow(X)
    q <- ncol(X)
    Interaction.mat <- X - matrix(1, nrow = p, ncol = 1) %*% 
        (matrix(1, nrow = 1, ncol = p) %*% X/p) - X %*% (rep(1, 
        q)/q) %*% matrix(1, nrow = 1, ncol = q) + mean(X)
    XMinOverallMean <- X - Xmean
    main.effects.rows <- apply(X, 1, mean) - mean(X)
    main.effects.cols <- apply(X, 2, mean) - mean(X)
    Residual.mat <- X - matrix(apply(X, 1, mean), ncol = 1) %*% 
        matrix(1, nrow = 1, ncol = ncol(X)) - matrix(1, nrow = nrow(X), 
        ncol = 1) %*% matrix(apply(X, 2, mean), nrow = 1) + mean(X)
    if (length(axis.col) < q) 
        axis.col <- (rep(axis.col, q))[1:q]
    if (length(ax.name.col) < q) 
        ax.name.col <- (rep(ax.name.col, q))[1:q]
    if (length(marker.col) < q) 
        marker.col <- (rep(marker.col, q))[1:q]
    if (length(tick.marker.col) < q) 
        tick.marker.col <- (rep(tick.marker.col, q))[1:q]
    svd.X <- svd(X)
    svd.XMinMean <- svd(XMinOverallMean)
    svd.InteractMat <- svd(Residual.mat)
    if (biad.variant == "Xmat") {
        Umat <- svd.X$u
        Sigma <- diag(svd.X$d)
        Vmat <- svd.X$v
        SVDmat <- X
        if (!is.null(X.new.rows)) 
            X.new.rows <- X.new.rows
        if (!is.null(X.new.columns)) 
            X.new.columns <- X.new.columns
    }
    if (biad.variant == "XminMeanMat") {
        Umat <- svd.XMinMean$u
        Sigma <- diag(svd.XMinMean$d)
        Vmat <- svd.XMinMean$v
        SVDmat <- XMinOverallMean
        if (!is.null(X.new.rows)) 
            X.new.rows <- X.new.rows - matrix(1, nrow = nrow(X.new.rows), 
                ncol = ncol(X.new.rows)) * Xmean
        if (!is.null(X.new.columns)) 
            X.new.columns <- X.new.columns - matrix(1, nrow = nrow(X.new.columns), 
                ncol = ncol(X.new.columns)) * Xmean
    }
    if (biad.variant == "InteractionMat") {
        Umat <- svd.InteractMat$u
        Sigma <- diag(svd.InteractMat$d)
        Vmat <- svd.InteractMat$v
        SVDmat <- Residual.mat
        if (!is.null(X.new.rows)) {
            X.new.rows.row.means <- apply(X.new.rows, 1, mean)
            X.column.means <- apply(X, 2, mean)
            X.new.rows <- X.new.rows - matrix(X.new.rows.row.means, 
                byrow = FALSE, nrow = nrow(X.new.rows), ncol = ncol(X)) - 
                matrix(X.column.means, byrow = TRUE, nrow = nrow(X.new.rows), 
                  ncol = ncol(X)) + matrix(1, nrow = nrow(X.new.rows), 
                ncol = ncol(X)) * Xmean
        }
        if (!is.null(X.new.columns)) {
            X.new.columns.column.means <- apply(X.new.columns, 
                2, mean)
            X.row.means <- apply(X, 1, mean)
            X.new.columns <- X.new.columns - matrix(X.row.means, 
                byrow = FALSE, nrow = nrow(X), ncol = ncol(X.new.columns)) - 
                matrix(X.new.columns.column.means, byrow = TRUE, 
                  nrow = nrow(X), ncol = ncol(X.new.columns)) + 
                matrix(1, nrow = nrow(X), ncol = ncol(X.new.columns)) * 
                  Xmean
        }
    }
    sing.values <- diag(Sigma)
    Quality <- sum((sing.values^2)[e.vects][1:dim.biplot])/sum((sing.values^2)) * 
        100
    Sigma2 <- Sigma %*% Sigma
    NumMat <- diag(diag(Vmat %*% Sigma2 %*% t(Vmat)))
    Uit.list.cols <- vector("list", q)
    Uit.list.Xhat <- vector("list", min(p, q))
    for (i in 1:min(p, q)) {
        Jr <- matrix(0, nrow = min(p, q), ncol = min(p, q))
        Jr[1:i, 1:i] <- diag(i)
        Uit.list.cols[[i]] <- diag(diag(Vmat %*% Sigma2 %*% Jr %*% 
            t(Vmat))) %*% solve(NumMat)
        Uit.list.Xhat[[i]] <- Umat %*% Sigma %*% Jr %*% t(Vmat)
        dimnames(Uit.list.Xhat[[i]]) <- dimnames(X)
    }
    Axis.predictivities <- round(diag(Uit.list.cols[[1]]), digits = 4)
    for (i in 2:min(p, q)) Axis.predictivities <- rbind(Axis.predictivities, 
        round(diag(Uit.list.cols[[i]]), digits = 4))
    dimnames(Axis.predictivities) <- list(paste("Dim", 1:min(p, 
        q), sep = "_"), dimnames(X)[[2]])
    if (predictivity.print) {
        if (dim.biplot == "1") 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], "(", 
                round(Axis.predictivities[1, ], digits = 2) * 
                  100, ")", sep = "")
        if (dim.biplot == "2") 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], "(", 
                round(Axis.predictivities[2, ], digits = 2) * 
                  100, ")", sep = "")
        if (dim.biplot == "3") 
            dimnames(X)[[2]] <- paste(dimnames(X)[[2]], "(", 
                round(Axis.predictivities[3, ], digits = 2) * 
                  100, ")", sep = "")
    }
    dev.new()
    options(pty = "s")
    par(mar = parplotmar)
    {
        vect.scaffolding <- e.vects[1:dim.biplot]
        if (SigmaHalf) {
            SigmaHalf <- Sigma^0.5
            SigmaMinHalf <- ifelse(SigmaHalf > 1e-09, 1/SigmaHalf, 
                0)
            XV <- (SVDmat %*% Vmat %*% SigmaMinHalf)[, vect.scaffolding] %*% 
                rotate.mat %*% reflect.mat
            V <- Vmat %*% SigmaHalf[, vect.scaffolding] %*% rotate.mat %*% 
                reflect.mat
            if (!is.null(X.new.rows)) 
                X.new.rows.plot <- (X.new.rows %*% Vmat %*% SigmaMinHalf)[, 
                  vect.scaffolding, drop = FALSE] %*% rotate.mat %*% 
                  reflect.mat
            if (!is.null(X.new.columns)) 
                X.new.columns.plot <- (t(X.new.columns) %*% Umat %*% 
                  SigmaMinHalf)[, vect.scaffolding, drop = FALSE] %*% 
                  rotate.mat %*% reflect.mat
        }
        else {
            XV <- (SVDmat %*% Vmat)[, vect.scaffolding] %*% rotate.mat %*% 
                reflect.mat
            V <- Vmat[, vect.scaffolding] %*% rotate.mat %*% 
                reflect.mat
            if (!is.null(X.new.rows)) 
                X.new.rows.plot <- (X.new.rows %*% Vmat)[, vect.scaffolding, 
                  drop = FALSE] %*% rotate.mat %*% reflect.mat
            if (!is.null(X.new.columns)) 
                X.new.columns.plot <- (t(X.new.columns) %*% Umat %*% 
                  Sigma^(-1))[, vect.scaffolding, drop = FALSE] %*% 
                  rotate.mat %*% reflect.mat
        }
        if (lambda) {
            lam.4 <- p * sum(V * V)/(q * sum(XV * XV))
            lam <- sqrt(sqrt(lam.4))
            XV <- XV * lam
            V <- V/lam
            if (!is.null(X.new.rows)) 
                X.new.rows.plot <- X.new.rows.plot * lam
            if (!is.null(X.new.columns)) 
                X.new.columns.plot <- X.new.columns.plot/lam
        }
        if (dim.biplot == 1) {
            XV <- cbind(XV, 0)
            V <- cbind(V, 0)
        }
        lims <- c(min(XV, V), max(XV, V))
        plot(rbind(XV[, 1:2], V[, 1:2]), asp = 1, xlim = lims * 
            exp.factor, ylim = lims * exp.factor, type = "n", 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        usr <- par("usr")
        plotshift <- propshift * (usr[4] - usr[3])
        if (dim.biplot == 1) {
            row.plot.coords <- cbind(XV[, 1, drop = FALSE], plotshift)
            dimnames(row.plot.coords)[[1]] <- dimnames(SVDmat)[[1]]
            col.plot.coords <- cbind(V[, 1, drop = FALSE], plotshift)
            if (!is.null(X.new.rows)) 
                X.new.rows.plot <- cbind(X.new.rows.plot[, 1, 
                  drop = FALSE], plotshift)
            if (!is.null(X.new.columns)) 
                X.new.columns.plot <- cbind(X.new.columns.plot[, 
                  1, drop = FALSE], plotshift)
        }
        else {
            row.plot.coords <- XV[, 1:2]
            col.plot.coords <- V
            if (!is.null(X.new.rows)) 
                X.new.rows.plot <- X.new.rows.plot[, 1:2, drop = FALSE]
            if (!is.null(X.new.columns)) 
                X.new.columns.plot <- X.new.columns.plot[, 1:2, 
                  drop = FALSE]
        }
        samples.coords <- XV
        points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        text(x = row.plot.coords, labels = dimnames(X)[[1]], 
            pos = text.pos[1], cex = row.points.label.size, col = row.points.col)
        if (!is.null(X.new.rows)) 
            points(X.new.rows.plot)
        if (!is.null(X.new.columns)) 
            points(X.new.columns.plot)
        if (plot.col.points) {
            points(col.plot.coords, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size, )
            if (column.points.text) 
                text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                  pos = text.pos[2], cex = column.points.label.size, 
                  col = column.points.col)
        }
        if (dim.biplot == 3) 
            dev.off()
        calibrations.list <- vector("list", q)
        ax.names <- dimnames(X)[[2]]
        maineffects.pos <- matrix(0, nrow = q, ncol = 2)
        for (i in 1:q) {
            calibrations.mat <- matrix(0)
            counter <- 0
            while (nrow(calibrations.mat) < 2 & counter < 11) {
                counter <- counter + 1
                markers.vals <- pretty(range(SVDmat[, i]), n = n.int[i])
                if (!is.null(expand.markers)) {
                  term <- markers.vals[1] - markers.vals[2]
                  for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                    markers.vals[length(markers.vals)] - term)
                  for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                    markers.vals[1] + term, after = 0)
                }
                transformed.markers.vals <- markers.vals
                if (add.maineffect & biad.variant == "InteractionMat") {
                  markers.vals <- pretty(range(SVDmat[, i]) + 
                    main.effects.cols[i], n = n.int[i])
                  if (!is.null(expand.markers)) {
                    term <- markers.vals[1] - markers.vals[2]
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[length(markers.vals)] - term)
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[1] + term, after = 0)
                  }
                  transformed.markers.vals <- markers.vals - 
                    main.effects.cols[i]
                }
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 1]
                  calibrations.y <- (transformed.markers.vals/sum(V[i, 
                    1:2]^2)) * V[i, 2]
                  if (show.maineffect) {
                    maineffect.x <- mean(SVDmat[, i])/sum(V[i, 
                      1:2]^2) * V[i, 1]
                    maineffect.y <- mean(SVDmat[, i])/sum(V[i, 
                      1:2]^2) * V[i, 2]
                    if (dim.biplot == 1) 
                      maineffect.y <- maineffect.y + i * constant + 
                        propshift
                    maineffects.pos[i, ] <- c(maineffect.x, maineffect.y)
                  }
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
            if (nrow(calibrations.mat) < 2) 
                warning(paste("Not 2 markers found on axis ", 
                  i, ". Increase n.int \n"))
            if (dim.biplot == 1 | dim.biplot == 2) {
                calibrations.x.in <- calibrations.mat[, 1]
                calibrations.y.in <- calibrations.mat[, 2]
                markers.vals.in <- calibrations.mat[, 3]
                calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                  calibrations.y.in, markers.vals.in), ncol = 3)
            }
            else calibrations.list[[i]] <- calibrations.mat
        }
        if (dim.biplot == 1 | dim.biplot == 2) {
            out <- (drawbipl.biad(calibrations.list = calibrations.list, 
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
            if (show.origin) 
                points(0, 0, pch = 3, cex = 1)
            if (!is.null(circle.proj)) {
                for (i in circle.proj) {
                  coordin <- row.plot.coords[i, 1:2]/2
                  rr <- sqrt(sum(coordin^2))
                  draw.circle(r = rr, h1 = coordin[1], h2 = coordin[2], 
                    col = "black")
                }
                if (show.maineffect) 
                  points(maineffects.pos, pch = 16, cex = 1)
            }
        }
        if (!select.origin & !is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval <- zoom(zoomval)
            interrupt <- FALSE
            plot(rbind(XV[, 1:2], V[, 1:2]), asp = 1, xlim = zoomval[1:2], 
                ylim = zoomval[3:4], type = "n", xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
                main = Title)
            usr <- par("usr")
            plotshift <- propshift * (usr[4] - usr[3])
            if (dim.biplot == 1) {
                row.plot.coords <- cbind(XV[, 1, drop = FALSE], 
                  plotshift)
                dimnames(row.plot.coords)[[1]] <- dimnames(SVDmat)[[1]]
                col.plot.coords <- cbind(V[, 1, drop = FALSE], 
                  plotshift)
                if (!is.null(X.new.rows)) 
                  X.new.rows.plot <- cbind(X.new.rows.plot[, 
                    1, drop = FALSE], plotshift)
                if (!is.null(X.new.columns)) 
                  X.new.columns.plot <- cbind(X.new.columns.plot[, 
                    1, drop = FALSE], plotshift)
            }
            else {
                row.plot.coords <- XV[, 1:2]
                col.plot.coords <- V
                if (!is.null(X.new.rows)) 
                  X.new.rows.plot <- X.new.rows.plot[, 1:2, drop = FALSE]
                if (!is.null(X.new.columns)) 
                  X.new.columns.plot <- X.new.columns.plot[, 
                    1:2, drop = FALSE]
            }
            samples.coords <- XV
            points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            text(x = row.plot.coords, labels = dimnames(X)[[1]], 
                pos = text.pos[1], cex = row.points.label.size, 
                col = row.points.col)
            if (!is.null(X.new.rows)) 
                points(X.new.rows.plot)
            if (!is.null(X.new.columns)) 
                points(X.new.columns.plot)
            if (plot.col.points) {
                points(col.plot.coords, pch = pch.col.points, 
                  col = column.points.col, cex = column.points.size, 
                  )
                if (column.points.text) 
                  text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                    pos = text.pos[2], cex = column.points.label.size, 
                    col = column.points.col)
            }
            if (dim.biplot == 3) 
                dev.off()
            calibrations.list <- vector("list", q)
            ax.names <- dimnames(X)[[2]]
            maineffects.pos <- matrix(0, nrow = q, ncol = 2)
            for (i in 1:q) {
                calibrations.mat <- matrix(0)
                counter <- 0
                while (nrow(calibrations.mat) < 2 & counter < 
                  11) {
                  counter <- counter + 1
                  markers.vals <- pretty(range(SVDmat[, i]), 
                    n = n.int[i])
                  if (!is.null(expand.markers)) {
                    term <- markers.vals[1] - markers.vals[2]
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[length(markers.vals)] - term)
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[1] + term, after = 0)
                  }
                  transformed.markers.vals <- markers.vals
                  if (add.maineffect & biad.variant == "InteractionMat") {
                    markers.vals <- pretty(range(SVDmat[, i]) + 
                      main.effects.cols[i], n = n.int[i])
                    if (!is.null(expand.markers)) {
                      term <- markers.vals[1] - markers.vals[2]
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[length(markers.vals)] - 
                          term)
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[1] + term, after = 0)
                    }
                    transformed.markers.vals <- markers.vals - 
                      main.effects.cols[i]
                  }
                  if (dim.biplot == 1 | dim.biplot == 2) {
                    calibrations.x <- (transformed.markers.vals/sum(V[i, 
                      1:2]^2)) * V[i, 1]
                    calibrations.y <- (transformed.markers.vals/sum(V[i, 
                      1:2]^2)) * V[i, 2]
                    if (show.maineffect) {
                      maineffect.x <- mean(SVDmat[, i])/sum(V[i, 
                        1:2]^2) * V[i, 1]
                      maineffect.y <- mean(SVDmat[, i])/sum(V[i, 
                        1:2]^2) * V[i, 2]
                      if (dim.biplot == 1) 
                        maineffect.y <- maineffect.y + i * constant + 
                          propshift
                      maineffects.pos[i, ] <- c(maineffect.x, 
                        maineffect.y)
                    }
                    calibrations.mat <- matrix(c(calibrations.x, 
                      calibrations.y, markers.vals), ncol = 3)
                    criterion.x <- calibrations.x > usr[1] & 
                      calibrations.x < usr[2]
                    criterion.y <- calibrations.y > usr[3] & 
                      calibrations.y < usr[4]
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
                    calibrations.mat <- cbind(calibrations.x, 
                      calibrations.y, calibrations.z, markers.vals)
                  }
                  n.int[i] <- 2 * n.int[i]
                }
                if (nrow(calibrations.mat) < 2) {
                  warning(paste("Not 2 markers found on axis ", 
                    i, ". Increase n.int \n"))
                  interrupt <- TRUE
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
            if (dim.biplot == 1 | dim.biplot == 2) {
                if (!interrupt) 
                  out <- (drawbipl.biad(calibrations.list = calibrations.list, 
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
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
                if (!is.null(circle.proj)) {
                  for (i in circle.proj) {
                    coordin <- row.plot.coords[i, 1:2]/2
                    rr <- sqrt(sum(coordin^2))
                    draw.circle(r = rr, h1 = coordin[1], h2 = coordin[2], 
                      col = "black")
                  }
                  if (show.maineffect) 
                    points(maineffects.pos, pch = 16, cex = 1)
                }
            }
        }
        if (select.origin) {
            cat("Move cursor where you want origin and press left mouse button \n")
            flush.console()
            origin.pos <- locator(1)
            calibrations.list <- vector("list", q)
            orthog.transx <- rep(origin.pos$x, q)
            orthog.transy <- rep(origin.pos$y, q)
            dev.new()
            n.int <- n.int.original
            options(pty = "s")
            par(mar = parplotmar)
            lims <- c(min(XV, V), max(XV, V))
            plot(rbind(XV[, 1:2], V[, 1:2]), asp = 1, xlim = lims * 
                exp.factor, ylim = lims * exp.factor, type = "n", 
                xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
                xaxs = "i", yaxs = "i", main = Title)
            usr <- par("usr")
            plotshift <- propshift * (usr[4] - usr[3])
            if (dim.biplot == 1) {
                row.plot.coords <- cbind(XV[, 1, drop = FALSE], 
                  plotshift)
                dimnames(row.plot.coords)[[1]] <- dimnames(SVDmat)[[1]]
                col.plot.coords <- cbind(V[, 1, drop = FALSE], 
                  plotshift)
                if (!is.null(X.new.rows)) 
                  X.new.rows.plot <- cbind(X.new.rows.plot[, 
                    1, drop = FALSE], plotshift)
                if (!is.null(X.new.columns)) 
                  X.new.columns.plot <- cbind(X.new.columns.plot[, 
                    1, drop = FALSE], plotshift)
            }
            else {
                row.plot.coords <- XV[, 1:2]
                col.plot.coords <- V
                if (!is.null(X.new.rows)) 
                  X.new.rows.plot <- X.new.rows.plot[, 1:2, drop = FALSE]
                if (!is.null(X.new.columns)) 
                  X.new.columns.plot <- X.new.columns.plot[, 
                    1:2, drop = FALSE]
            }
            samples.coords <- XV
            points(row.plot.coords, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            text(x = row.plot.coords, labels = dimnames(X)[[1]], 
                pos = text.pos[1], cex = row.points.label.size, 
                col = row.points.col)
            if (!is.null(X.new.rows)) 
                points(X.new.rows.plot)
            if (!is.null(X.new.columns)) 
                points(X.new.columns.plot)
            col.plot.coords <- V
            if (plot.col.points) {
                points(col.plot.coords, pch = pch.col.points, 
                  col = column.points.col, cex = column.points.size, 
                  )
                if (column.points.text) 
                  text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                    pos = text.pos[2], cex = column.points.label.size, 
                    col = column.points.col)
            }
            points(0, 0, pch = 3, cex = 1)
            maineffects.pos <- matrix(0, nrow = q, ncol = 2)
            for (i in 1:q) {
                calibrations.mat <- matrix(0)
                phi.vec <- diag(1/diag(V %*% t(V))) %*% V %*% 
                  c(orthog.transx[i], orthog.transy[i])
                counter <- 0
                while (nrow(calibrations.mat) < 2 & counter < 
                  9) {
                  counter <- counter + 1
                  markers.vals <- pretty(range(SVDmat[, i]), 
                    n = n.int[i])
                  if (!is.null(expand.markers)) {
                    term <- markers.vals[1] - markers.vals[2]
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[length(markers.vals)] - term)
                    for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                      markers.vals[1] + term, after = 0)
                  }
                  transformed.markers.vals <- markers.vals
                  if (add.maineffect & biad.variant == "InteractionMat") {
                    markers.vals <- pretty(range(SVDmat[, i]) + 
                      main.effects.cols[i], n = n.int[i])
                    if (!is.null(expand.markers)) {
                      term <- markers.vals[1] - markers.vals[2]
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[length(markers.vals)] - 
                          term)
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[1] + term, after = 0)
                    }
                    transformed.markers.vals <- markers.vals - 
                      main.effects.cols[i]
                  }
                  if (dim.biplot == 1 | dim.biplot == 2) {
                    calibrations.x <- (transformed.markers.vals/sum(V[i, 
                      1:2]^2)) * V[i, 1]
                    calibrations.x <- orthog.transx[i] + calibrations.x - 
                      phi.vec[i] * V[i, 1]
                    calibrations.y <- (transformed.markers.vals/sum(V[i, 
                      1:2]^2)) * V[i, 2]
                    calibrations.y <- orthog.transy[i] + calibrations.y - 
                      phi.vec[i] * V[i, 2]
                    if (show.maineffect) {
                      maineffect.x <- mean(SVDmat[, i])/sum(V[i, 
                        1:2]^2) * V[i, 1]
                      maineffect.x <- orthog.transx[i] + maineffect.x - 
                        phi.vec[i] * V[i, 1]
                      maineffect.y <- mean(SVDmat[, i])/sum(V[i, 
                        1:2]^2) * V[i, 2]
                      maineffect.y <- orthog.transy[i] + maineffect.y - 
                        phi.vec[i] * V[i, 2]
                      maineffects.pos[i, ] <- c(maineffect.x, 
                        maineffect.y)
                    }
                    calibrations.mat <- matrix(c(calibrations.x, 
                      calibrations.y, markers.vals), ncol = 3)
                    criterion.x <- calibrations.x > usr[1] & 
                      calibrations.x < usr[2]
                    criterion.y <- calibrations.y > usr[3] & 
                      calibrations.y < usr[4]
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
                    calibrations.mat <- cbind(calibrations.x, 
                      calibrations.y, calibrations.z, markers.vals)
                  }
                  n.int[i] <- 2 * n.int[i]
                }
                if (nrow(calibrations.mat) < 2) 
                  warning(paste("Not 2 markers found on axis ", 
                    i, ". Increase n.int \n"))
                if (dim.biplot == 1 | dim.biplot == 2) {
                  calibrations.x.in <- calibrations.mat[, 1]
                  calibrations.y.in <- calibrations.mat[, 2]
                  markers.vals.in <- calibrations.mat[, 3]
                  calibrations.list[[i]] <- matrix(c(calibrations.x.in, 
                    calibrations.y.in, markers.vals.in), ncol = 3)
                }
                else calibrations.list[[i]] <- calibrations.mat
            }
            {
                out <- (drawbipl.biad(calibrations.list = calibrations.list, 
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
                if (show.origin) 
                  points(0, 0, pch = 3, cex = 1)
                if (!is.null(circle.proj)) {
                  for (j in circle.proj) {
                    new.center.x <- min(row.plot.coords[j, 1], 
                      orthog.transx[1]) + abs(row.plot.coords[j, 
                      1] - orthog.transx[1])/2
                    new.center.y <- min(row.plot.coords[j, 2], 
                      orthog.transy[1]) + abs(row.plot.coords[j, 
                      2] - orthog.transy[1])/2
                    rr <- (sqrt(((row.plot.coords[j, 1] - orthog.transx[1])^2) + 
                      ((row.plot.coords[j, 2] - orthog.transy[1])^2)))/2
                    draw.circle(r = rr, h1 = new.center.x, h2 = new.center.y, 
                      col = "black")
                  }
                }
                if (show.maineffect) 
                  points(maineffects.pos, pch = 16, cex = 1)
            }
            if (!is.null(zoomval)) {
                if (!is.numeric(zoomval)) 
                  stop("zoomval must be numeric")
                cat("Indicate zooming position \n")
                flush.console()
                zoomval <- zoom(zoomval)
                interrupt <- FALSE
                plot(rbind(XV[, 1:2], V[, 1:2]), asp = 1, xlim = zoomval[1:2], 
                  ylim = zoomval[3:4], type = "n", xlab = "", 
                  ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                  yaxs = "i", main = Title)
                usr <- par("usr")
                plotshift <- propshift * (usr[4] - usr[3])
                if (dim.biplot == 1) {
                  row.plot.coords <- cbind(XV[, 1, drop = FALSE], 
                    plotshift)
                  dimnames(row.plot.coords)[[1]] <- dimnames(SVDmat)[[1]]
                  col.plot.coords <- cbind(V[, 1, drop = FALSE], 
                    plotshift)
                  if (!is.null(X.new.rows)) 
                    X.new.rows.plot <- cbind(X.new.rows.plot[, 
                      1, drop = FALSE], plotshift)
                  if (!is.null(X.new.columns)) 
                    X.new.columns.plot <- cbind(X.new.columns.plot[, 
                      1, drop = FALSE], plotshift)
                }
                else {
                  row.plot.coords <- XV[, 1:2]
                  col.plot.coords <- V
                  if (!is.null(X.new.rows)) 
                    X.new.rows.plot <- X.new.rows.plot[, 1:2, 
                      drop = FALSE]
                  if (!is.null(X.new.columns)) 
                    X.new.columns.plot <- X.new.columns.plot[, 
                      1:2, drop = FALSE]
                }
                samples.coords <- XV
                points(row.plot.coords, pch = pch.row.points, 
                  col = row.points.col, cex = row.points.size)
                text(x = row.plot.coords, labels = dimnames(X)[[1]], 
                  pos = text.pos[1], cex = row.points.label.size, 
                  col = row.points.col)
                if (!is.null(X.new.rows)) 
                  points(X.new.rows.plot)
                if (!is.null(X.new.columns)) 
                  points(X.new.columns.plot)
                col.plot.coords <- V
                if (plot.col.points) {
                  points(col.plot.coords, pch = pch.col.points, 
                    col = column.points.col, cex = column.points.size, 
                    )
                  if (column.points.text) 
                    text(x = col.plot.coords, labels = dimnames(X)[[2]], 
                      pos = text.pos[2], cex = column.points.label.size, 
                      col = column.points.col)
                }
                points(0, 0, pch = 3, cex = 1)
                maineffects.pos <- matrix(0, nrow = q, ncol = 2)
                for (i in 1:q) {
                  calibrations.mat <- matrix(0)
                  phi.vec <- diag(1/diag(V %*% t(V))) %*% V %*% 
                    c(orthog.transx[i], orthog.transy[i])
                  counter <- 0
                  while (nrow(calibrations.mat) < 2 & counter < 
                    9) {
                    counter <- counter + 1
                    markers.vals <- pretty(range(SVDmat[, i]), 
                      n = n.int[i])
                    if (!is.null(expand.markers)) {
                      term <- markers.vals[1] - markers.vals[2]
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[length(markers.vals)] - 
                          term)
                      for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                        markers.vals[1] + term, after = 0)
                    }
                    transformed.markers.vals <- markers.vals
                    if (add.maineffect & biad.variant == "InteractionMat") {
                      markers.vals <- pretty(range(SVDmat[, i]) + 
                        main.effects.cols[i], n = n.int[i])
                      if (!is.null(expand.markers)) {
                        term <- markers.vals[1] - markers.vals[2]
                        for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                          markers.vals[length(markers.vals)] - 
                            term)
                        for (jj in 1:expand.markers) markers.vals <- append(markers.vals, 
                          markers.vals[1] + term, after = 0)
                      }
                      transformed.markers.vals <- markers.vals - 
                        main.effects.cols[i]
                    }
                    if (dim.biplot == 1 | dim.biplot == 2) {
                      calibrations.x <- (transformed.markers.vals/sum(V[i, 
                        1:2]^2)) * V[i, 1]
                      calibrations.x <- orthog.transx[i] + calibrations.x - 
                        phi.vec[i] * V[i, 1]
                      calibrations.y <- (transformed.markers.vals/sum(V[i, 
                        1:2]^2)) * V[i, 2]
                      calibrations.y <- orthog.transy[i] + calibrations.y - 
                        phi.vec[i] * V[i, 2]
                      if (show.maineffect) {
                        maineffect.x <- mean(SVDmat[, i])/sum(V[i, 
                          1:2]^2) * V[i, 1]
                        maineffect.x <- orthog.transx[i] + maineffect.x - 
                          phi.vec[i] * V[i, 1]
                        maineffect.y <- mean(SVDmat[, i])/sum(V[i, 
                          1:2]^2) * V[i, 2]
                        maineffect.y <- orthog.transy[i] + maineffect.y - 
                          phi.vec[i] * V[i, 2]
                        maineffects.pos[i, ] <- c(maineffect.x, 
                          maineffect.y)
                      }
                      calibrations.mat <- matrix(c(calibrations.x, 
                        calibrations.y, markers.vals), ncol = 3)
                      criterion.x <- calibrations.x > usr[1] & 
                        calibrations.x < usr[2]
                      criterion.y <- calibrations.y > usr[3] & 
                        calibrations.y < usr[4]
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
                      calibrations.mat <- cbind(calibrations.x, 
                        calibrations.y, calibrations.z, markers.vals)
                    }
                    n.int[i] <- 2 * n.int[i]
                  }
                  if (nrow(calibrations.mat) < 2) {
                    warning(paste("Not 2 markers found on axis ", 
                      i, ". Increase n.int \n"))
                    interrupt <- TRUE
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
                {
                  if (!interrupt) 
                    out <- (drawbipl.biad(calibrations.list = calibrations.list, 
                      row.plot.coords.mat = row.plot.coords, 
                      ax = ax, axes.names = ax.names, ax.name.size = ax.name.size, 
                      axis.col = axis.col, ax.name.col = ax.name.col, 
                      constant = constant, line.length = line.length, 
                      markers = markers, marker.size = marker.size, 
                      marker.col = marker.col, offset = offset, 
                      offset.m = offset.m, ort.lty = ort.lty, 
                      pos = pos, pos.m = pos.m, propshift = propshift, 
                      predictions.sample = predictions.sample, 
                      predictions.allsamples.onaxis = predictions.allsamples.onaxis, 
                      side.label = side.label, tick.marker.col = tick.marker.col, 
                      q = q))
                  if (show.origin) 
                    points(0, 0, pch = 3, cex = 1)
                  if (!is.null(circle.proj)) {
                    for (j in circle.proj) {
                      new.center.x <- min(row.plot.coords[j, 
                        1], orthog.transx[1]) + abs(row.plot.coords[j, 
                        1] - orthog.transx[1])/2
                      new.center.y <- min(row.plot.coords[j, 
                        2], orthog.transy[1]) + abs(row.plot.coords[j, 
                        2] - orthog.transy[1])/2
                      rr <- (sqrt(((row.plot.coords[j, 1] - orthog.transx[1])^2) + 
                        ((row.plot.coords[j, 2] - orthog.transy[1])^2)))/2
                      draw.circle(r = rr, h1 = new.center.x, 
                        h2 = new.center.y, col = "black")
                    }
                  }
                  if (show.maineffect) 
                    points(maineffects.pos, pch = 16, cex = 1)
                }
            }
        }
        if (dim.biplot == 3) {
            dimnames(col.plot.coords)[[1]] <- dimnames(SVDmat)[[2]]
            out <- (drawbipl.3dim.biad(Zrows = samples.coords, 
                Zcols = col.plot.coords, z.axes = calibrations.list, 
                X.new.rows = points.new, X.new.rows.pars = X.new.rows.pars, 
                axes.names = ax.names, adj.3d = adj.3d, alpha.3d = alpha.3d, 
                ax.col = axis.col, ax.col.3d = ax.col.3d, aspect.3d = aspect.3d, 
                cex.3d = cex.3d, col.plane.3d = col.plane.3d, 
                column.points.col = column.points.col, row.points.col = row.points.col, 
                plot.col.points = plot.col.points, col.text.3d = col.text.3d, 
                factor.x = factor.x, factor.y = factor.y, font.3d = font.3d, 
                ID.labs = ID.labs, ID.3d = ID.3d, samples.plot = samples.plot, 
                size.ax.3d = size.ax.3d, size.points.3d = size.points.3d, 
                Titles.3d = Titles.3d))
            if (predictions.3D) {
                predictions <- biad.predictivities(X = X, e.vects = e.vects, 
                  biad.variant = biad.variant, add.maineffect = add.maineffect, 
                  predictions.dim = 3)$predictions.rows
                dimnames(predictions) <- dimnames(X)
            }
            else predictions <- NULL
        }
    }
    if (!lambda) 
        lam <- NULL
    if (legend.show) {
        dev.new()
        par(mar = rep(0.25, 4))
        plot(1:10, 1:10, typ = "n", xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", bty = "n")
        usr <- par("usr")
        legend(x = usr[1], y = usr[4], legend = dimnames(X)[[1]], 
            pch = pch.row.points, col = row.points.col, ncol = legend.columns[1], 
            text.col = row.points.col)
        legend(x = usr[1], y = usr[4]/2, legend = dimnames(X)[[2]], 
            pch = pch.col.points, col = column.points.col, ncol = legend.columns[2], 
            lty = 1, text.col = column.points.col)
    }
    out.list <- list(out, X.mat = X, lambda = lam, Quality = Quality, 
        predictions = predictions, Axis.predictivities = Axis.predictivities, 
        main.effects.cols = main.effects.cols, main.effects.rows = main.effects.rows, 
        SVDmat = SVDmat, Umat = Umat, Vmat = Vmat, Sigma = Sigma)
    out.list[output]
}
