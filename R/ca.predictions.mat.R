ca.predictions.mat <-
function (Pmat, e.vects = 1:2) 
{
    r.mass <- apply(Pmat, 1, sum)
    c.mass <- apply(Pmat, 2, sum)
    Dr <- diag(r.mass)
    Dc <- diag(c.mass)
    R.profile.mat <- solve(Dr) %*% Pmat
    dimnames(R.profile.mat) <- dimnames(Pmat)
    Smat <- sqrt(solve(Dr)) %*% (Pmat - matrix(r.mass, ncol = 1) %*% 
        matrix(c.mass, nrow = 1)) %*% sqrt(solve(Dc))
    svd.Smat <- svd(Smat)
    Phi.mat <- sqrt(solve(Dr)) %*% svd.Smat$u
    GH.mat <- sqrt(solve(Dc)) %*% svd.Smat$v
    F.mat <- Phi.mat %*% diag(svd.Smat$d)
    row.profile.predictions <- (F.mat[, e.vects, drop = FALSE] %*% 
        t(GH.mat[, e.vects, drop = FALSE]) + 1) %*% Dc
    dimnames(row.profile.predictions) <- dimnames(Pmat)
    row.profile.predictions
}
