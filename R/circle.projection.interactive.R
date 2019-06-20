circle.projection.interactive <-
function (origin = c(0, 0), colr = "black", cent.line = FALSE, 
    ...) 
{
    loc <- locator(1)
    coordin <- c(loc$x, loc$y)/2
    origin <- origin/2
    rr <- sqrt(sum((coordin - origin)^2))
    draw.circle(r = rr, h1 = coordin[1] + origin[1], h2 = coordin[2] + 
        origin[2], col = colr, ...)
    if (cent.line) 
        lines(x = c(origin[1], coordin[1]) * 2, y = c(origin[2], 
            coordin[2]) * 2, col = colr)
}
