MonoPlot.skew <-
function (X = Rothkopf.vowels, form.triangle1 = NULL, form.triangle2 = NULL, 
    ...) 
{
    if (nrow(X) != ncol(X)) 
        stop("X must be a square matrix")
    K.vec <- rep(c(1, -1), nrow(X)%/%2)
    if (!identical(nrow(X)%%2, 0)) 
        K.vec <- c(K.vec, 1)
    K <- diag(K.vec)
    M <- 0.5 * (X + t(X))
    N <- 0.5 * (X - t(X))
    svd.uit <- svd(N)
    u.vec <- svd.uit$u[, 1, drop = FALSE]
    v.vec <- svd.uit$u[, 2, drop = FALSE]
    par(pty = "s")
    plot(u.vec, v.vec, asp = 1, pty = "n", xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    points(u.vec, v.vec, pch = 16, col = "black", ...)
    text(u.vec, v.vec, labels = dimnames(X)[[1]], pos = c(1, 
        1, 1, 2, 1), cex = 0.9)
    if (!is.null(form.triangle1)) 
        polygon(c(0, u.vec[form.triangle1]), c(0, v.vec[form.triangle1]), 
            col = "green")
    if (!is.null(form.triangle1)) 
        polygon(c(0, u.vec[form.triangle2]), c(0, v.vec[form.triangle2]), 
            col = "red")
    points(0, 0, pch = 3, cex = 5, col = "darkgrey")
    sigma <- diag(svd.uit$d)
    list(M = M, N = N, K = K, V = svd.uit$v, U = svd.uit$u, sigma = diag(svd.uit$d), 
        KtU = K %*% t(svd.uit$u))
}
