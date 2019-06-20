cabipl.doubling <-
function (X, N, ...) 
{
    p <- nrow(X)
    q <- ncol(X)
    mat.double <- cbind(X, N * matrix(1, nrow = p, ncol = q) - 
        X)
    colnames(mat.double) <- c(paste(colnames(X), "plus", sep = "."), 
        paste(colnames(X), "min", sep = "."))
    cabipl(mat.double, ...)
}
