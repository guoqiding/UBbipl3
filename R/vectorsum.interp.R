vectorsum.interp <-
function (vertex.points = 3, p = vertex.points, from = c(0, 0), 
    pch.centroid = 15, col.centroid = "black", size.centroid = 1, 
    pch.interp = 16, col.interp = "red", col.arrow = "red", size.interp = 1.2, 
    length = 0.15, angle = 15, lty = 1, lwd = 1.75, arrow.to.centroid.only = FALSE, 
    point = FALSE, ...) 
{
    out <- locator(vertex.points)
    polygon(x = out$x, y = out$y, ...)
    n <- vertex.points
    xvec <- out$x
    yvec <- out$y
    xvec.0 <- c(out$x[n], out$x[-n])
    yvec.0 <- c(out$y[n], out$y[-n])
    centroid.x <- mean(xvec)
    centroid.y <- mean(yvec)
    points(x = centroid.x, y = centroid.y, pch = pch.centroid, 
        col = col.centroid, cex = size.centroid)
    end.x <- p * (centroid.x - from[1]) + from[1]
    end.y <- p * (centroid.y - from[2]) + from[2]
    arrows(from[1], from[2], centroid.x, centroid.y, length = length, 
        angle = angle, lty = 1, col = col.centroid, lwd = lwd + 
            0.25)
    if (point) 
        points(x = centroid.x, y = centroid.y, pch = pch.centroid, 
            col = col.centroid, cex = size.centroid)
    if (!arrow.to.centroid.only) {
        arrows(centroid.x, centroid.y, end.x, end.y, length = length, 
            angle = angle, lty = lty, col = col.arrow, lwd = lwd)
        if (point) 
            points(x = end.x, y = end.y, pch = pch.interp, col = col.interp, 
                cex = size.interp)
    }
    list(out = out, centroid.x = centroid.x, centroid.y = centroid.y)
}
