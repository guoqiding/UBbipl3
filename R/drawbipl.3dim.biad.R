drawbipl.3dim.biad <-
function (Zrows, Zcols, z.axes, axes.names, X.new.rows, X.new.rows.pars, 
    adj.3d = 0.5, alpha.3d = 0.7, ax.col = rep("black", ncol(Z)), 
    ax.col.3d = "black", aspect.3d = "iso", cex.3d = 0.6, col.plane.3d = "lightgrey", 
    col.text.3d = "black", factor.x = 2, factor.y = 2, font.3d = 2, 
    ID.labs = FALSE, ID.3d = dimnames(Zrows)[[1]], samples.plot = TRUE, 
    plot.col.points = TRUE, size.ax.3d = 1, size.points.3d = 50, 
    row.points.col, column.points.col, Titles.3d = c("", "", 
        "Dim 1", "Dim 2", "Dim 3")) 
{
    require(rgl)
    z.axes.names = axes.names
    graph.coord <- z.axes[[1]][c(1, nrow(z.axes[[1]])), 1:3]
    for (i in 2:length(z.axes)) graph.coord <- rbind(graph.coord, 
        z.axes[[i]][c(1, nrow(z.axes[[i]])), 1:3])
    open3d()
    texts3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
    points3d(graph.coord, size = 0)
    aspect3d("iso")
    if (samples.plot) 
        for (i in 1:nrow(Zrows)) spheres3d(Zrows[i, 1:3, drop = FALSE], 
            size = size.points.3d, color = row.points.col[i], 
            alpha = alpha.3d)
    if (ID.labs) 
        texts3d(Zrows[, 1:3], text = ID.3d, adj = adj.3d, col = row.points.col, 
            font = font.3d, cex = cex.3d)
    if (plot.col.points) 
        for (i in 1:nrow(Zcols)) points3d(Zcols[i, 1:3, drop = FALSE], 
            size = size.points.3d, color = column.points.col[i])
    for (i in 1:length(z.axes)) {
        axis <- z.axes[[i]]
        texts3d(axis[axis[, 4] == max(axis[, 4]), 1:3], text = z.axes.names[i], 
            col = col.text.3d[i], font = font.3d, cex = cex.3d)
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
    device.ID <- rgl.cur()
    answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    while (!(answer == "Y" | answer == "y" | answer == "N" | 
        answer == "n")) answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    if (answer == "Y" | answer == "y") 
        repeat {
            file.name <- readline("Provide file name including full NOT in quotes and SINGLE back slashes!  \n")
            file.name <- paste(file.name, ".png", sep = "")
            snapshot3d(file = file.name)
            rgl.set(device.ID)
            answer2 <- readline("Save another 3D graph as a .png file? Y/N  \n")
            if (answer2 == "Y" | answer2 == "y") 
                next
            else break
        }
    else rgl.set(device.ID)
}
