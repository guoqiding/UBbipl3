MDSbipl <-
function (X = as.matrix(Ocotea.data[, 3:8]), scaled.mat = FALSE, 
    iter = 500, dist.metric = "Pythagoras", rep.init = 2, axis.col = rep("darkgrey", 
        ncol(X)), ax.name.col = rep("black", ncol(X)), ax.name.size = 0.8, 
    calibrate = TRUE, marker.size = 0.6, marker.col = rep("darkgrey", 
        ncol(X)), n.int = rep(5, ncol(X)), offset = c(0, 0.2, 
        0.1, 0.1), offset.m = rep(0, ncol(X)), pch.samples = 15, 
    pch.samples.col = c(rep("red", 20), rep("blue", 7), rep("green", 
        10)), pch.samples.size = 1.2, pos = "Hor", pos.m = rep(2, 
        ncol(X)), side.label = rep("left", ncol(X)), strepie = c(1, 
        1)) 
{
    n <- nrow(X)
    p <- ncol(X)
    X.original <- X
    if (scaled.mat) 
        X <- scale(X)
    else X <- scale(X, center = TRUE, scale = FALSE)
    if (dist.metric == "Pythagoras") 
        Dmat <- sqrt(Pythagoras.dist(X))
    W.mat <- matrix(0, nrow = n, ncol = n)
    W.mat[lower.tri(W.mat, diag = T)] <- 1
    W.mat <- t(W.mat)
    Ini.mat <- matrix(rnorm(2 * n), nrow = n)
    out.select <- SMACOF(Dmat, W.mat, Ini.mat, iter = iter)
    sigma.r.select <- out.select$sigma.r
    for (i in 1:rep.init) {
        Ini.mat <- matrix(rnorm(2 * n), nrow = n)
        out <- SMACOF(Dmat, W.mat, Ini.mat, iter = iter, scaled.mat = scaled.mat)
        if (out$sigma.r < sigma.r.select) {
            out.select <- out
            sigma.r.select <- out$sigma.r
        }
    }
    Zmat <- out.select$X.up
    Bmat <- solve(t(Zmat) %*% Zmat) %*% t(Zmat) %*% X
    rownames(Zmat) <- rownames(X)
    par(pty = "s")
    plot(Zmat, type = "p", asp = 1, col = pch.samples.col, pch = pch.samples, 
        cex = pch.samples.size, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "")
    text(Zmat, labels = rownames(Zmat), pos = 1, cex = 0.6, offset = 0.2)
    marker.list <- vector("list", p)
    axes.rows <- solve(diag(diag(t(Bmat) %*% Bmat))) %*% t(Bmat)
    for (i in 1:p) marker.list[[i]] <- find.markers(jj = i, X = X.original, 
        axes.rows = axes.rows, n.int = n.int[i])
    if (calibrate) {
        usr <- par("usr")
        for (i in 1:p) {
            mark.mat <- marker.list[[i]]
            marker.mat.invals <- mark.mat[(mark.mat[, 1] >= usr[1] & 
                mark.mat[, 1] <= usr[2] & mark.mat[, 2] >= usr[3] & 
                mark.mat[, 2] <= usr[4]), ]
            gradient <- Draw.line2(x.vals = marker.mat.invals[, 
                1], y.vals = marker.mat.invals[, 2], marker.vals = marker.mat.invals[, 
                3], line.name = colnames(X)[i], offset = offset, 
                pos = pos, axis.col = axis.col, ax.name.col = ax.name.col, 
                ax.name.size = ax.name.size)$gradient
            marker.mat.invals <- matrix(marker.mat.invals[marker.mat.invals[, 
                4] == 1, ], ncol = 4)
            for (j in 1:length(marker.mat.invals[, 1])) {
                Draw.onecmline(x = marker.mat.invals[j, 1], y = marker.mat.invals[j, 
                  2], grad = -1/gradient, expand = strepie[1], 
                  both.sides = TRUE, col = marker.col[i])
                Plot.marker(x = marker.mat.invals[j, 1], y = marker.mat.invals[j, 
                  2], grad = -1/gradient, mark = marker.mat.invals[j, 
                  3], expand = strepie[2], marker.size = marker.size, 
                  col = marker.col[i], pos.m = pos.m[i], offset.m = offset.m[i], 
                  side.label = side.label[i])
            }
        }
    }
    list(iterasies = out.select$iterasies, Bmat = Bmat, axes.rows = axes.rows, 
        sigma.r = sigma.r.select, sum.Dmat = sum(Dmat), X.up = out.select$X.up)
}
