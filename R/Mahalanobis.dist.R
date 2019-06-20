Mahalanobis.dist <-
function (X, S.inv = solve(var(X))) 
{
    swd <- svd(S.inv)
    S.inv.sqrt <- swd$u %*% diag(swd$d^0.5) %*% t(swd$u)
    Pythagoras.dist(X %*% S.inv.sqrt)
}
