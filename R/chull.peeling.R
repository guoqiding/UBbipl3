chull.peeling <-
function (datmat, k = 1, colours = UBcolours, line.width = 1.5, 
    ...) 
{
    plot(datmat[, 1:2], type = "p", ...)
    for (i in (1:k)) {
        points <- chull(datmat)
        lines(datmat[c(points, points[1]), ], col = colours[i], 
            lwd = line.width)
        datmat <- datmat[-points, ]
    }
}
