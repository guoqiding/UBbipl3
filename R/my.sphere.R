my.sphere <-
function (datmat, width, minvec, colour = "red") 
{
    datmat <- scale(datmat, center = minvec, scale = FALSE)
    spheres3d(datmat[, 1], datmat[, 3], datmat[, 2], radius = 10/width, 
        col = colour)
}
