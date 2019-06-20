drawbipl.ca <-
function (row.plot.coords.mat, calibrations.list, ax, axes.names = NULL, 
    ax.name.size = 0.5, axis.col = "grey", ax.name.col = "grey", 
    constant, line.length = c(1, 1), markers = TRUE, marker.size = 0.5, 
    marker.col = "grey", offset, offset.m, ort.lty, pos, pos.m, 
    predictions.sample = NULL, propshift, side.label, tick.marker.col = "grey", 
    q, X.new.points = NULL) 
{
    for (i in 1:length(calibrations.list)) {
        if (!is.matrix(calibrations.list[[i]])) 
            stop("calibration.mat must be a matrix \n")
        if (ncol(calibrations.list[[i]]) != 3) 
            stop("calibration.mat must have three columns \n")
        if (nrow(calibrations.list[[i]]) < 2) 
            stop("calibration.mat must at least two rows. Increase n.int \n")
    }
    if (is.null(ax)) 
        axes <- NULL
    else axes <- (1:q)[ax]
    if (!is.null(predictions.sample)) {
        predictions <- data.frame(matrix(NA, nrow = length(axes), 
            ncol = length(predictions.sample)))
        dimnames(predictions) <- list(1:nrow(predictions), paste("s", 
            predictions.sample, sep = ""))
    }
    if (!is.null(X.new.points)) {
        predictions.new.points <- matrix(NA, nrow = length(axes), 
            ncol = nrow(X.new.points))
        dimnames(predictions.new.points) <- list(1:nrow(predictions.new.points), 
            paste("np", 1:ncol(predictions.new.points), sep = "."))
    }
    else predictions.new.points <- NULL
    r.names <- NULL
    for (i in axes) {
        calibrations.x.in <- calibrations.list[[i]][, 1]
        calibrations.y.in <- calibrations.list[[i]][, 2] + i * 
            constant + propshift
        markers.vals.in <- calibrations.list[[i]][, 3]
        if (is.null(axes.names)) 
            axis.name <- paste("col", i, sep = "")
        else axis.name <- axes.names[i]
        r.names[i] <- axis.name
        axis.col.1 <- axis.col[i]
        ax.name.col.1 <- ax.name.col[i]
        out <- Draw.line2(x.vals = calibrations.x.in, y.vals = calibrations.y.in, 
            marker.vals = markers.vals.in, line.name = axis.name, 
            ax.name.size = ax.name.size, axis.col = axis.col.1, 
            ax.name.col = ax.name.col.1, offset = offset, pos = pos)
        if (!is.null(predictions.sample)) {
            for (ss in predictions.sample) {
                predictions[i, paste("s", ss, sep = "")] <- round(DrawOrthogline(x1 = calibrations.x.in[1], 
                  y1 = calibrations.y.in[1], x2 = calibrations.x.in[length(calibrations.x.in)], 
                  y2 = calibrations.y.in[length(calibrations.y.in)], 
                  val1 = markers.vals.in[1], val2 = markers.vals.in[length(calibrations.x.in)], 
                  px = row.plot.coords.mat[ss, 1], py = row.plot.coords.mat[ss, 
                    2], ort.lty = ort.lty), digits = 6)
            }
        }
        if (!is.null(X.new.points)) {
            for (ss in 1:nrow(X.new.points)) {
                predictions.new.points[i, ss] <- round(DrawOrthogline(x1 = calibrations.x.in[1], 
                  y1 = calibrations.y.in[1], x2 = calibrations.x.in[length(calibrations.x.in)], 
                  y2 = calibrations.y.in[length(calibrations.y.in)], 
                  val1 = markers.vals.in[1], val2 = markers.vals.in[length(calibrations.x.in)], 
                  px = X.new.points[ss, 1], py = X.new.points[ss, 
                    2], ort.lty = ort.lty), digits = 6)
            }
        }
        gradient <- out$gradient
        for (j in 1:length(calibrations.x.in)) Draw.onecmline(x = calibrations.x.in[j], 
            y = calibrations.y.in[j], grad = -1/gradient, expand = line.length[1], 
            both.sides = TRUE, col = tick.marker.col)
        if (markers == TRUE) {
            x.labvals <- calibrations.x.in
            y.labvals <- calibrations.y.in
            tick.labels <- zapsmall(markers.vals.in)
        }
        else {
            x.labvals <- calibrations.x.in[c(1, length(calibrations.x.in))]
            y.labvals <- calibrations.y.in[c(1, length(calibrations.y.in))]
            tick.labels <- zapsmall(markers.vals.in)
            tick.labels <- tick.labels[c(1, length(tick.labels))]
        }
        for (j in 1:length(x.labvals)) Plot.marker(x = x.labvals[j], 
            y = y.labvals[j], grad = -1/gradient, mark = tick.labels[j], 
            expand = line.length[2], marker.size = marker.size, 
            col = marker.col, pos.m = pos.m[i], offset.m = offset.m[i], 
            side.label = side.label[i])
    }
    if (is.null(predictions.sample)) 
        predictions <- NULL
    if (!is.null(predictions)) {
        predictions <- na.omit(predictions)
        dimnames(predictions)[[1]] <- r.names[!is.na(r.names)]
    }
    if (!is.null(predictions.new.points)) {
        predictions.new.points <- na.omit(predictions.new.points)
        dimnames(predictions.new.points)[[1]] <- r.names[!is.na(r.names)]
    }
    if (is.null(predictions.new.points)) 
        return(predictions)
    else return(predictions.new.points)
}
