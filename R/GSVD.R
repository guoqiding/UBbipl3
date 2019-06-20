GSVD <-
function (A, Omega, Fi, Kster = 2) 
{
    svd.O <- svd(Omega)
    svd.F <- svd(Fi)
    Omega.vkw <- svd.O$u %*% sqrt(diag(svd.O$d)) %*% t(svd.O$u)
    Fi.vkw <- svd.F$u %*% sqrt(diag(svd.F$d)) %*% t(svd.F$u)
    BB <- (Omega.vkw) %*% A %*% Fi.vkw
    swd <- svd(BB)
    U.mat <- solve(Omega.vkw) %*% swd$u[, 1:Kster]
    V.mat <- solve(Fi.vkw) %*% swd$v[, 1:Kster]
    Lambda.mat <- diag(swd$d[1:Kster])
    list(Kster = Kster, U = U.mat, V = V.mat, Lambda = Lambda.mat)
}
