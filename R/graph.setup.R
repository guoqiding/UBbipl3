graph.setup <-
function (x.axis = c(0, 1), y.axis = c(0, 1), z.axis = c(0, 1)) 
{
    open3d()
    aspect3d("iso")
    bg3d("#FFFFFF", fogtype = "lin")
    view3d(theta = 200, phi = 25, fov = 1)
    if (x.axis[2] < x.axis[1]) 
        x.axis <- rev(x.axis)
    if (y.axis[2] < y.axis[1]) 
        y.axis <- rev(y.axis)
    if (z.axis[2] < z.axis[1]) 
        z.axis <- rev(z.axis)
    x.min <- x.axis[1]
    y.min <- y.axis[1]
    z.min <- z.axis[1]
    width <- max(x.axis[2] - x.axis[1], y.axis[2] - y.axis[1], 
        z.axis[2] - z.axis[1])
    if (x.axis[2] - x.axis[1] < width) 
        x.axis[2] <- x.axis[1] + width
    if (y.axis[2] - y.axis[1] < width) 
        y.axis[2] <- y.axis[1] + width
    if (z.axis[2] - z.axis[1] < width) 
        z.axis[2] <- z.axis[1] + width
    points3d(x.axis - x.min, z.axis - z.min, y.axis - y.min, 
        alpha = 0)
    bbox3d(col = "#FEFEFE", alpha = 0.25, xat = 0, yat = 0, zat = 0)
    list(width, c(x.min, y.min, z.min))
}
