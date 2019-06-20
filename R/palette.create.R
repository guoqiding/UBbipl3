palette.create <-
function (mat = NULL) 
{
    if (is.null(mat)) 
        mat <- matrix(c(228, 55, 77, 152, 255, 255, 166, 247, 
            153, 26, 126, 175, 78, 127, 255, 86, 129, 153, 28, 
            184, 74, 163, 0, 51, 40, 191, 153), ncol = 3)
    rgb(mat[, 1], mat[, 2], mat[, 3], max = 255)
}
