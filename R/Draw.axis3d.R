Draw.axis3d <-
function (direction, markers, axis, O.vec, minvec, width = 0.05, 
    marker.width = c(0.8, 0.8, 0.8), colour = "black", bothsides = TRUE, 
    begin.pos = markers[1], mark.labels = rep(TRUE, length(markers)), 
    lines.size = 2, markers.size = 1) 
{
    axis.no <- (1:3)[c("x", "y", "z") == axis]
    axis.mean <- O.vec[axis.no]
    mark.pos <- markers - axis.mean
    if (bothsides) 
        axis.begin <- direction * (min(mark.pos) - (mark.pos[2] - 
            mark.pos[1])/1.1)
    else axis.begin <- direction * (begin.pos - axis.mean)
    axis.end <- direction * (max(mark.pos) + (mark.pos[2] - mark.pos[1])/1.1)
    axis.begin <- axis.begin + O.vec - minvec
    axis.end <- axis.end + O.vec - minvec
    my.arrow(c(axis.begin[1], axis.end[1]), c(axis.begin[3], 
        axis.end[3]), c(axis.begin[2], axis.end[2]), size.g = "large", 
        colour = colour, size = 3, width = width)
    if (bothsides) 
        my.arrow(c(axis.end[1], axis.begin[1]), c(axis.end[3], 
            axis.begin[3]), c(axis.end[2], axis.begin[2]), size.g = "large", 
            colour = colour, size = 3, width = width)
    marker.basis <- get.orthog(direction, width)
    for (m in 1:length(markers)) {
        marker.pos <- direction * mark.pos[m]
        axis.markers <- marker.pos + O.vec - minvec
        my.marker <- marker.basis + cbind(axis.markers, axis.markers)
        lines3d(my.marker[1, ], my.marker[3, ], my.marker[2, 
            ], col = colour, size = lines.size)
        axis.markers <- axis.markers - marker.width
        axis.markers <- my.marker[, 2] - marker.width
        if (mark.labels[m]) 
            text3d(axis.markers[1], axis.markers[3], axis.markers[2], 
                markers[m], col = colour, cex = markers.size)
    }
    direction %*% t(mark.pos)
}
