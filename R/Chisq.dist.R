Chisq.dist <-
function (X) 
{
    R.mat <- diag(apply(X, 1, sum))
    C.mat <- diag(apply(X, 2, sum))
    RMinEenXCMinHalf <- solve(R.mat) %*% X %*% solve(sqrt(C.mat))
    B.mat <- RMinEenXCMinHalf %*% t(RMinEenXCMinHalf)
    Sq.Chis2.dis <- matrix(1, nrow = nrow(X), ncol = nrow(X)) %*% 
        diag(diag(B.mat)) + t(matrix(1, nrow = nrow(X), ncol = nrow(X)) %*% 
        diag(diag(B.mat))) - 2 * B.mat
    Sq.Euclid.dis <- matrix(1, nrow = nrow(X), ncol = nrow(X)) %*% 
        diag(diag(X %*% t(X))) + t(matrix(1, nrow = nrow(X), 
        ncol = nrow(X)) %*% diag(diag(X %*% t(X)))) - 2 * X %*% 
        t(X)
    list(Sq.Chis2.dis = Sq.Chis2.dis, Sq.Euclid.dis = Sq.Euclid.dis)
}
