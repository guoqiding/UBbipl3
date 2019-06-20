my.points <-
function (datmat, width, minvec) 
{
    datmat <- scale(datmat, center = minvec, scale = FALSE)
    points3d(datmat[, 1], datmat[, 3], datmat[, 2], size = 150/width, 
        col = colour.fun(0, 0, 255))
}
