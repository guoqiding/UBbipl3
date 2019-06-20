draw.polygon <-
function (vertex.points = 3, p = vertex.points, ...) 
{
    out <- locator(vertex.points)
    polygon(x = out$x, y = out$y, ...)
    n <- vertex.points
    xvec <- out$x
    yvec <- out$y
    xvec.0 <- c(out$x[n], out$x[-n])
    yvec.0 <- c(out$y[n], out$y[-n])
    Area <- (sum(xvec.0 * yvec - xvec * yvec.0))/2
    centroid.x <- (sum((xvec + xvec.0) * (xvec.0 * yvec - xvec * 
        yvec.0)))/(6 * Area)
    centroid.y <- (sum((yvec + yvec.0) * (xvec.0 * yvec - xvec * 
        yvec.0)))/(6 * Area)
    list(out = out, centroid.x = centroid.x, centroid.y = centroid.y)
}
