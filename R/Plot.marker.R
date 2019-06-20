Plot.marker <-
function (x, y, grad, mark, expand = 10, marker.size, col = 1, 
    offset.m = 0.5, pos.m = 1, side.label = "left") 
{
    usr <- par("usr")
    uin <- par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
    mm <- 1/(uin[1] * 25.4)
    d <- (expand + 2) * mm
    b <- d * sqrt(1/(1 + grad * grad))
    a <- b * grad
    if (is.infinite(grad)) {
        if (sign(grad) > 0) 
            text(x = x, y = y + d, labels = mark, cex = marker.size, 
                pos = pos.m, offset = offset.m, col = col)
        else text(x = x, y = y - d, labels = mark, cex = marker.size, 
            pos = pos.m, offset = offset.m, col = col)
    }
    else {
        if (side.label == "right") 
            text(x = x + b, y = y + a, labels = mark, cex = marker.size, 
                pos = pos.m, offset = offset.m, col = col)
        if (side.label == "left") 
            text(x = x - b, y = y - a, labels = mark, cex = marker.size, 
                pos = pos.m, offset = offset.m, col = col)
    }
}
