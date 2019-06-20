my.star <-
function (xpos = 0, ypos = 0, size = 1, colour = "red", points = 5) 
{
    degrees <- (0:(points * 2 - 1)) * (360/(points * 2))
    angles <- degrees * pi/180
    radii <- rep(c(0.5, 1), points) * size
    x.vec <- radii * cos(angles)
    y.vec <- radii * sin(angles)
    polygon(x.vec + xpos, y.vec + ypos, border = NA, col = colour)
}
