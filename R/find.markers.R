find.markers <-
function (jj, axes.rows, X, n.int) 
{
    means <- apply(X, 2, mean)
    sds <- sqrt(apply(X, 2, var))
    number.points <- 20
    std.markers <- pretty(X[, jj], n = n.int)
    std.range <- c(min(std.markers), max(std.markers))
    std.markers.min <- std.markers - (std.range[2] - std.range[1])
    std.markers.max <- std.markers + (std.range[2] - std.range[1])
    std.markers <- c(std.markers, std.markers.min, std.markers.max)
    interval <- (std.markers - means[jj])/sds[jj]
    axis.vals <- seq(from = min(interval), to = max(interval), 
        length = number.points)
    axis.vals <- sort(unique(c(axis.vals, interval)))
    number.points <- length(axis.vals)
    axis.points <- matrix(0, nrow = number.points, ncol = 4)
    axis.points[, 1] <- (axis.vals) * axes.rows[jj, 1]
    axis.points[, 2] <- (axis.vals) * axes.rows[jj, 2]
    axis.points[, 3] <- axis.vals * sds[jj] + means[jj]
    axis.points[, 4] <- 0
    for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
        3] - std.markers) == 0)) 
        axis.points[i, 4] <- 1
    axis.points
}
