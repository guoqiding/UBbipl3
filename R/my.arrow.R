my.arrow <-
function (x, y, z, width = 0.1, colour = "black", size.g = "small", 
    ...) 
{
    lines3d(x, y, z, col = colour, ...)
    vector <- c(x[2] - x[1], y[2] - y[1], z[2] - z[1])
    up.point <- c(x[2], y[2], z[2]) + 0.01 * vector
    length <- sqrt(sum(vector^2))
    width <- min(width, length/50)
    if (size.g == "small") 
        basis.middle <- up.point - 7 * width * vector
    else basis.middle <- up.point - 0.85 * width * vector
    basis <- get.orthog(vector, width)
    basis <- basis + cbind(basis.middle, basis.middle)
    triangles3d(c(up.point[1], basis[1, ]), c(up.point[2], basis[2, 
        ]), c(up.point[3], basis[3, ]), col = colour, alpha = 1, 
        front = "fill", back = "fill", size = 2)
}
