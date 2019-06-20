draw.circle.line <-
function (r = 1, h1 = 0, h2 = 0, len = 1000, theta = NULL, calib = NULL, 
    pos.m = 2, offset.m = 0, ...) 
{
    a <- seq(from = 0, to = 2 * pi, len = len)
    Y <- cbind(r * cos(a) + h1, r * sin(a) + h2)
    lines(Y, ...)
    if (!is.null(theta)) 
        if (!is.null(calib)) {
            xcoords <- seq(from = h1, to = r * cos(theta) + h1, 
                len = calib + 1)
            ycoords <- seq(from = h2, to = r * sin(theta) + h2, 
                len = calib + 1)
            intcept <- ycoords + xcoords/tan(theta)
            labs <- c("", round((1:calib)/calib, digits = 2))
            for (i in 1:length(intcept)) {
                Draw.onecmline.label(x = xcoords[i], y = ycoords[i], 
                  grad = -1/tan(theta), both.sides = TRUE, expand = 1, 
                  label = NULL)
                Plot.marker(x = xcoords[i], y = ycoords[i], grad = -1/tan(theta), 
                  expand = 1, mark = labs[i], marker.size = 0.5, 
                  col = 1, offset = offset.m, pos = pos.m)
            }
        }
}
