draw.arrow <-
function (arrow = 2, ...) 
{
    out <- locator(arrow)
    arrows(x0 = out$x[1], x1 = out$x[2], y0 = out$y[1], y1 = out$y[2], 
        ...)
}
