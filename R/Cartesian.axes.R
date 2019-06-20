Cartesian.axes <-
function (x.axis, y.axis, z.axis, minvec, width) 
{
    x.markers <- pretty(x.axis)
    x.markers <- x.markers[x.markers > minvec[1]]
    y.markers <- pretty(y.axis)
    y.markers <- y.markers[y.markers > minvec[2]]
    z.markers <- pretty(z.axis)
    z.markers <- z.markers[z.markers > minvec[3]]
    Draw.axis3d(c(1, 0, 0), x.markers, axis = "x", O.vec = minvec, 
        minvec = minvec, bothsides = F, begin.pos = x.axis[1])
    text3d(x.axis[2] + 45/width, 0, -100/width, "X", col = "black")
    Draw.axis3d(c(0, 1, 0), y.markers, axis = "y", O.vec = minvec, 
        minvec = minvec, bothsides = F, begin.pos = y.axis[1])
    text3d(-70/width, 0, y.axis[2] + 45/width, "Y", col = "black")
    Draw.axis3d(c(0, 0, 1), z.markers, axis = "z", O.vec = minvec, 
        minvec = minvec, bothsides = F, begin.pos = z.axis[1], 
        marker.width = c(0.8 - 0.6, 0.8 - 3.5, 0.8))
    text3d(0, z.axis[2] + 30/width, 0, "Z", col = "black")
}
