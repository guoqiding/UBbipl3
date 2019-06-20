plane.intersect <-
function (Vmat, mu, x.val, y.axis, minvec, line.col = colour.fun(255, 
    204, 0), line.width) 
{
    xy <- expand.grid(x.val - mu[1], y.axis - mu[2])
    z.vec <- Vmat[3, 2]/Vmat[1, 2] * xy[, 1] + ((Vmat[3, 2] * 
        Vmat[1, 1] - Vmat[3, 1] * Vmat[1, 2])/(Vmat[2, 2] * Vmat[1, 
        1] - Vmat[2, 1] * Vmat[1, 2])) * (xy[, 2] - xy[, 1] * 
        Vmat[2, 2]/Vmat[1, 2])
    lines3d(rep(x.val - minvec[1], 2), z.vec + mu[3] - minvec[3], 
        y.axis - minvec[2], col = line.col, lwd = line.width)
}
