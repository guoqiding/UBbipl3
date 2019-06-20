draw.circle <-
function (r = 1, h1 = 0, h2 = 0, len = 1000, ...) 
{
    a <- seq(from = 0, to = 2 * pi, len = len)
    Y <- cbind(r * cos(a) + h1, r * sin(a) + h2)
    lines(Y, ...)
}
