MonoPlot.coefvar <-
function (X = as.matrix(Flotation.data[, -1]), as.axes = FALSE, 
    axis.col = "darkgrey", ax.name.col = "black", ax.name.size = 0.65, 
    calibrate = FALSE, dim.plot = 2, lambda = FALSE, marker.size = 0.75, 
    marker.col = rep(1, ncol(X)), n.int = rep(5, ncol(X)), offset = c(0.5, 
        0.5, 0.5, 0.5), offset.m = rep(0, ncol(X)), pos = "Orthog", 
    pos.m = rep(1, ncol(X)), samples.plot = FALSE, side.label = rep("left", 
        ncol(X)), strepie = c(1, 1), VJ.plot = FALSE, zoomval = NULL) 
{
    n <- nrow(X)
    p <- ncol(X)
    means <- apply(X, 2, mean)
    delta.mat <- diag(means)
    print(delta.mat)
    coefvar.vec <- apply(X, 2, sd)/apply(X, 2, mean)
    names(coefvar.vec) <- colnames(X)
    E.mat <- (diag(n) - matrix(1, nrow = n, ncol = n)/n) %*% 
        X %*% solve(delta.mat)
    svd.out <- svd(E.mat)
    lam <- 1
    Jmat <- diag(c(rep(1, dim.plot), rep(0, min(n, p) - dim.plot)))
    V <- svd.out$v
    U <- svd.out$u
    Sigma <- diag(svd.out$d)
    VSigma <- V %*% Sigma
    if (lambda) {
        lam.4 <- n * sum(VSigma * VSigma)/(p * sum(U * U))
        lam <- sqrt(sqrt(lam.4))
        U <- U * lam
        VSigma <- VSigma/lam
    }
    if (VJ.plot) 
        vars.points <- V %*% Jmat
    else vars.points <- VSigma %*% Jmat
    par(pty = "s")
    marker.mat <- vector("list", p)
    X.mat <- X
    for (ii in 1:p) marker.mat[[ii]] <- find.markers(jj = ii, 
        X = X.mat, axes.rows = vars.points[, 1:2, drop = FALSE], 
        n.int = n.int[ii])
    sample.points <- U %*% Sigma %*% Jmat
    if (samples.plot) 
        plot(rbind(vars.points, sample.points), type = "n", asp = 1, 
            xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    else plot(vars.points, type = "n", asp = 1, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    if (as.axes) {
        for (i in 1:p) {
            theta <- atan(coefficients(lm(c(0, vars.points[i, 
                2]) ~ c(0, vars.points[i, 1])))[2])
        }
        a <- par("usr")
        for (i in 1:p) {
            mark.mat <- marker.mat[[i]]
            marker.mat.invals <- mark.mat[(mark.mat[, 1] >= a[1] & 
                mark.mat[, 1] <= a[2] & mark.mat[, 2] >= a[3] & 
                mark.mat[, 2] <= a[4]), ]
            gradient <- Draw.line2(x.vals = marker.mat.invals[, 
                1], y.vals = marker.mat.invals[, 2], marker.vals = marker.mat.invals[, 
                3], line.name = dimnames(X)[[2]][i], offset = offset, 
                pos = pos, axis.col = axis.col, ax.name.col = ax.name.col, 
                ax.name.size = ax.name.size)$gradient
            marker.mat.invals <- matrix(marker.mat.invals[marker.mat.invals[, 
                4] == 1, ], ncol = 4)
            if (calibrate) {
                for (j in 1:length(marker.mat.invals[, 1])) Plot.marker(x = marker.mat.invals[j, 
                  1], y = marker.mat.invals[j, 2], grad = -1/gradient, 
                  mark = marker.mat.invals[j, 3], expand = strepie[2], 
                  marker.size = marker.size, col = marker.col[i], 
                  pos.m = pos.m[i], offset.m = offset.m[i], side.label = side.label[i])
            }
        }
    }
    points(vars.points, pch = 16, col = "red", asp = 1, cex = 1.5)
    if (samples.plot) 
        points(sample.points, pch = 15, col = "green")
    if (!is.null(zoomval)) {
        if (!is.numeric(zoomval)) 
            stop("zoomval must be numeric")
        zoomval <- zoom(zoomval)
        if (samples.plot) 
            plot(rbind(vars.points, sample.points), xlim = zoomval[1:2], 
                ylim = zoomval[3:4], type = "n", asp = 1, xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "")
        else plot(vars.points, xlim = zoomval[1:2], ylim = zoomval[3:4], 
            type = "n", asp = 1, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "")
        if (as.axes) {
            for (i in 1:p) {
                theta <- atan(coefficients(lm(c(0, vars.points[i, 
                  2]) ~ c(0, vars.points[i, 1])))[2])
            }
            a <- par("usr")
            for (i in 1:p) {
                mark.mat <- marker.mat[[i]]
                marker.mat.invals <- mark.mat[(mark.mat[, 1] >= 
                  a[1] & mark.mat[, 1] <= a[2] & mark.mat[, 2] >= 
                  a[3] & mark.mat[, 2] <= a[4]), ]
                gradient <- Draw.line2(x.vals = marker.mat.invals[, 
                  1], y.vals = marker.mat.invals[, 2], marker.vals = marker.mat.invals[, 
                  3], line.name = dimnames(X)[[2]][i], offset = offset, 
                  pos = pos, axis.col = axis.col, ax.name.col = ax.name.col, 
                  ax.name.size = ax.name.size)$gradient
                marker.mat.invals <- matrix(marker.mat.invals[marker.mat.invals[, 
                  4] == 1, ], ncol = 4)
                if (calibrate) {
                  for (j in 1:length(marker.mat.invals[, 1])) Plot.marker(x = marker.mat.invals[j, 
                    1], y = marker.mat.invals[j, 2], grad = -1/gradient, 
                    mark = marker.mat.invals[j, 3], expand = strepie[2], 
                    marker.size = marker.size, col = marker.col[i], 
                    pos.m = pos.m[i], offset.m = offset.m[i], 
                    side.label = side.label[i])
                }
            }
        }
        points(vars.points, pch = 16, col = "red", asp = 1, cex = 1.5)
        if (samples.plot) 
            points(sample.points, pch = 15, col = "green")
    }
    adequacies <- diag((svd.out$v[, 1:2]) %*% t(svd.out$v[, 1:2])/(svd.out$v) %*% 
        t(svd.out$v))
    names(adequacies) <- colnames(X)
    points(x = 0, y = 0, pch = 3, cex = 4)
    list(cov.X = cov(X), cor.X = cor(X), coefvar.vec = coefvar.vec, 
        adequacies = adequacies, svd.out$d^2, svd.out$v, svd.out$v %*% 
            diag(svd.out$d), VSigma)
}
