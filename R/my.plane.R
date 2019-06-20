my.plane <-
function (Vmat, mu, x.axis, z.axis, minvec, colour = c(colour.fun(255, 
    255, 0), colour.fun(95, 51, 0))) 
{
    xz <- expand.grid(x.axis, z.axis)
    y.vec <- Vmat[2, 2]/Vmat[1, 2] * xz[, 1] + ((Vmat[2, 2] * 
        Vmat[1, 1] - Vmat[1, 2] * Vmat[2, 1])/(Vmat[3, 2] * Vmat[1, 
        1] - Vmat[3, 1] * Vmat[1, 2])) * (xz[, 2] - xz[, 1] * 
        Vmat[3, 2]/Vmat[1, 2])
    surface3d(x.axis + mu[1] - minvec[1], z.axis + mu[3] - minvec[3], 
        matrix(y.vec, ncol = 2) + mu[2] - minvec[2], col = colour, 
        alpha = 0.5)
}
