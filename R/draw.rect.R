draw.rect <-
function (vertex.points = 2, ...) 
{
    out <- locator(vertex.points)
    rect(xleft = out$x[1], xright = out$x[2], ybottom = out$y[1], 
        ytop = out$y[2], ...)
}
