ConCentrEllipse <-
function (X, kappa = 2, covmat = NULL, col = 1, type = c("l", 
    "p", "b"), ...) 
{
    type <- type[1]
    means <- matrix(apply(X, 2, mean), nrow = 2)
    if (is.null(covmat)) 
        covmat <- cov(X)
    else covmat <- covmat
    range.vec <- apply(X, 2, range)
    mid.vec <- apply(range.vec, 2, function(x) (x[2] + x[1])/2)
    dif <- max(range.vec[2, ] - range.vec[1, ])/2
    xlim <- c(mid.vec[1] - dif, mid.vec[1] + dif)
    ylim <- c(mid.vec[2] - dif, mid.vec[2] + dif)
    if (type == "p" | type == "b") {
        par(pty = "s")
        plot(x = X[, 1], y = X[, 2], xlab = "x", ylab = "y", 
            xlim = xlim, ylim = ylim, col = col, asp = 1, ...)
    }
    svd.covmat <- svd(covmat)
    a <- (0:6283)/1000
    Y <- cbind(cos(a), sin(a))
    Y <- Y %*% diag(sqrt(svd.covmat$d)) %*% t(svd.covmat$v) * 
        kappa
    Y <- Y + matrix(rep(1, 6284), ncol = 1) %*% t(means)
    lines(Y, col = col, ...)
}
