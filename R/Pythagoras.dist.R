Pythagoras.dist <-
function (X) 
{
    n <- nrow(X)
    B <- X %*% t(X)
    Delta <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) for (j in 1:n) Delta[i, j] <- B[i, i] + B[j, 
        j] - 2 * B[i, j]
    Delta
}
