Eigen.twosided <-
function (A, B, eps = 1e-08) 
{
    if (!(all((A - t(A)) <= 1e-06))) 
        stop("A not symmetric")
    if (!(all((B - t(B)) <= 1e-06))) 
        stop("B not symmetric")
    if (any(eigen(B)$values < eps)) 
        stop("B not positive definite")
    svd.B.out <- svd(B)
    B.sqrt <- svd.B.out$u %*% diag(sqrt(svd.B.out$d)) %*% t(svd.B.out$u)
    svd.2.out <- svd(solve(B.sqrt) %*% A %*% solve(B.sqrt))
    W <- solve(B.sqrt) %*% svd.2.out$u
    Lambda.mat <- diag(svd.2.out$d)
    list(A.mat = A, B.mat = B, W.mat = W, Lambda.mat = Lambda.mat)
}
