MonoPlot.cor <-
function (X = Flotation.data[, -1], as.axes = TRUE, axis.col = "darkgrey", 
    ax.name.col = "black", ax.name.size = 0.65, arrows = TRUE, 
    calib.no = 5, circle = TRUE, dim.plot = 2, lambda = FALSE, 
    n.int = rep(5, ncol(X)), offset = c(0.5, 0.5, 0.5, 0.5), 
    offset.m = rep(0, ncol(X)), plot.vars.points = TRUE, pos = "Orthog", 
    pos.m = rep(1, ncol(X)), print.ax.approx = FALSE) 
{
    n <- nrow(X)
    p <- ncol(X)
    lam <- 1
    X.norm <- scale(X)/sqrt(n - 1)
    svd.out <- svd(X.norm)
    Jmat <- diag(c(rep(1, dim.plot), rep(0, min(n, p) - dim.plot)))
    U <- svd.out$u
    V.cor <- svd.out$v
    Sigma.cor <- diag(svd.out$d)
    VSigma.cor <- V.cor %*% Sigma.cor
    if (identical(lambda, TRUE)) {
        lam.4 <- n * sum(VSigma.cor * VSigma.cor)/(p * sum(U * 
            U))
        lam <- sqrt(sqrt(lam.4))
        U <- U * lam
        VSigma.cor <- VSigma.cor/lam
    }
    vars.points.cor <- VSigma.cor %*% Jmat
    par(pty = "s")
    marker.mat <- vector("list", p)
    X.mat <- X
    for (ii in 1:p) marker.mat[[ii]] <- find.markers(jj = ii, 
        X = X.mat, axes.rows = vars.points.cor[, 1:2, drop = FALSE], 
        n.int = n.int[ii])
    plot(VSigma.cor %*% Jmat, type = "n", asp = 1, xaxt = "n", 
        yaxt = "n", xlab = "", ylab = "")
    if (as.axes) {
        theta <- rep(NA, p)
        theta2 <- rep(NA, p)
        for (i in 1:p) {
            theta[i] <- atan(coefficients(lm(c(0, vars.points.cor[i, 
                2]) ~ c(0, vars.points.cor[i, 1])))[2])
            theta2[i] <- atan2(vars.points.cor[i, 2], vars.points.cor[i, 
                1])
        }
        if (circle) {
            for (i in 1:p) {
                draw.circle.line(r = 1/lam, theta = theta2[i], 
                  calib = calib.no, col = "green", pos.m = pos.m[i], 
                  offset.m = offset.m[i])
            }
        }
    }
    unit.cor.approx <- sqrt(diag(VSigma.cor %*% Jmat %*% t(VSigma.cor %*% 
        Jmat)))
    names(unit.cor.approx) <- colnames(X)
    a <- par("usr")
    if (print.ax.approx) 
        axis.name <- paste(colnames(X), " (", round(unit.cor.approx, 
            digits = 2), ")", sep = "")
    else axis.name <- colnames(X)
    for (i in 1:p) {
        mark.mat <- marker.mat[[i]]
        marker.mat.invals <- mark.mat[(mark.mat[, 1] >= a[1] & 
            mark.mat[, 1] <= a[2] & mark.mat[, 2] >= a[3] & mark.mat[, 
            2] <= a[4]), ]
        Draw.line2(x.vals = marker.mat.invals[, 1], y.vals = marker.mat.invals[, 
            2], marker.vals = marker.mat.invals[, 3], line.name = axis.name[i], 
            offset = offset, pos = pos, axis.col = axis.col, 
            ax.name.col = ax.name.col, ax.name.size = ax.name.size)
    }
    if (identical(arrows, TRUE)) 
        arrows(x0 = rep(0, p), y0 = rep(0, p), x1 = cos(theta2) * 
            unit.cor.approx, y1 = sin(theta2) * unit.cor.approx, 
            col = "red", lwd = 1.8)
    if (plot.vars.points) 
        points(vars.points.cor, col = "red", pch = 16, cex = 1.2)
    predictivities <- diag((svd.out$v %*% diag(svd.out$d)[, 1:2]) %*% 
        t(svd.out$v %*% diag(svd.out$d)[, 1:2])/(svd.out$v %*% 
        diag(svd.out$d)) %*% t(svd.out$v %*% diag(svd.out$d)))
    adequacies <- diag((svd.out$v[, 1:2]) %*% t(svd.out$v[, 1:2])/((svd.out$v) %*% 
        t(svd.out$v)))
    list(cor.X = cor(X), adequacies = adequacies, predictivities = predictivities, 
        unit.cor.approx = unit.cor.approx, VSigma.cor %*% Jmat %*% 
            t(VSigma.cor %*% Jmat), VSigma.cor %*% t(VSigma.cor))
}
