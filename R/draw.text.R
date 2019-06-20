draw.text <-
function (vertex.points = 1, string, xpd = NA, ...) 
{
    out <- locator(vertex.points)
    text(x = out$x[1], y = out$y[1], labels = string, xpd = xpd, 
        ...)
}
