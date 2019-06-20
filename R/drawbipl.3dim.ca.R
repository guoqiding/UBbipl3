drawbipl.3dim.ca <-
function (Z, z.axes, axes.names, X.new, X.new.pars, alpha.3d = 0.7, 
    ax.col = rep("black", ncol(Z)), adjust.3d = c(0.5, 0.5), 
    ax.col.3d = "black", aspect.3d = "iso", cex.3d = 0.6, col.samples = "green", 
    col.plane.3d = "lightgrey", col.text.3d = "black", factor.x = 2, 
    factor.y = 2, font.3d = 2, ID.labs = FALSE, ID.3d = 1:nrow(Z), 
    samples.plot = TRUE, size.ax.3d = 1, size.points.3d = 10, 
    Titles.3d = c("", "", "Dim 1", "Dim 2", "Dim 3")) 
{
    require(rgl)
    z.axes.names = axes.names
    graph.coord <- z.axes[[1]][c(1, nrow(z.axes[[1]])), 1:3]
    for (i in 2:length(z.axes)) graph.coord <- rbind(graph.coord, 
        z.axes[[i]][c(1, nrow(z.axes[[i]])), 1:3])
    open3d()
    texts3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
    points3d(graph.coord, size = 0)
    aspect3d(aspect.3d)
    if (samples.plot) 
        for (i in 1:nrow(Z)) spheres3d(Z[i, 1:3, drop = FALSE], 
            radius = size.points.3d, color = col.samples)
    if (ID.labs) 
        texts3d(Z[, 1:3], adj = adjust.3d, texts = ID.3d, col = col.samples, 
            font = font.3d, cex = cex.3d)
    for (i in 1:length(z.axes)) {
        axis <- z.axes[[i]]
        texts3d(axis[axis[, 4] == max(axis[, 4]), 1:3], text = z.axes.names[i], 
            col = col.text.3d, font = font.3d, cex = cex.3d)
        texts3d(axis[axis[, 4] == min(axis[, 4]), 1:3], text = z.axes.names[i], 
            col = col.text.3d, font = font.3d, cex = cex.3d)
        axis <- axis[axis[, 4] != max(axis[, 4]), 1:3]
        lines3d(axis[c(1, nrow(axis)), 1:3], size = size.ax.3d, 
            color = ax.col[i])
    }
    surface3d(x = c(min(graph.coord[, 1]), max(graph.coord[, 
        1])) * factor.x, y = c(min(graph.coord[, 2]), max(graph.coord[, 
        2])) * factor.y, z = matrix(0, nrow = 2, ncol = 2), col = col.plane.3d, 
        alpha = alpha.3d)
    axes3d(col = ax.col.3d)
    title3d(main = Titles.3d[1], sub = Titles.3d[2], xlab = Titles.3d[3], 
        ylab = Titles.3d[4], zlab = Titles.3d[5])
    if (!(is.null(X.new))) {
        rgl.points(X.new, size = X.new.pars[[3]], color = X.new.pars[[2]])
        rgl.texts(X.new, text = X.new.pars[[4]], adj = c(0.5, 
            0.5), col = X.new.pars[[2]], font = font.3d, cex = X.new.pars[[5]])
    }
    device.ID <- rgl.cur()
    answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    while (!(answer == "Y" | answer == "y" | answer == "N" | 
        answer == "n")) answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    if (answer == "Y" | answer == "y") 
        repeat {
            file.name <- readline("Provide file name including full NOT in quotes and SINGLE back slashes!  \n")
            file.name <- paste(file.name, ".png", sep = "")
            rgl.snapshot(file = file.name)
            rgl.set(device.ID)
            answer2 <- readline("Save another 3D graph as a .png file? Y/N  \n")
            if (answer2 == "Y" | answer2 == "y") 
                next
            else break
        }
    else rgl.set(device.ID)
}
