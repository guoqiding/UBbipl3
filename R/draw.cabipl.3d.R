draw.cabipl.3d <-
function (Z, z.axes, z.axes.names, X.new.profiles = NULL, e.vects = 1:3, 
    alpha.3d = 0.7, ax.col, ax.col.3d = "black", aspect.3d = "iso", 
    cex.3d = 0.6, col.plane.3d = "lightgrey", col.text.3d = "black", 
    class.vec, factor.x = 2, factor.y = 2, font.3d = 2, GH.mat, 
    ID.labs = FALSE, ID.3d = 1:nrow(Z), p, plot.class.means = TRUE, 
    size.ax.3d = 0.5, size.means.3d = 10, size.points.3d = 5, 
    specify.classes, Title = NULL, Titles.3d = c("", "", "x", 
        "y", "z"), Z.means.mat) 
{
    {
        require(rgl)
        graph.coord <- z.axes[[1]][c(1, nrow(z.axes[[1]])), 1:3]
        plot(1:10, 1:10, type = "n", xlab = "", ylab = "", yaxt = "n", 
            xaxt = "n", bty = "n")
        legend(x = 5, y = 5, legend = dimnames(Z.means.mat)[[1]], 
            pch = 15, col = Z.means.mat[, 4], cex = 1.5)
        .samples.plot.3d <- function(Z, Z.means.mat, class.vec, 
            specify.classes, size.points.3d) {
            if (is.null(Z.means.mat)) 
                J <- 0
            else J <- nrow(Z.means.mat)
            classes <- unique(dimnames(Z.means.mat)[[1]])[specify.classes]
            for (j in classes) {
                Z.class <- data.frame(Z[class.vec == j, , drop = FALSE], 
                  stringsAsFactors = FALSE)
                points3d(Z.class[, 1:3], size = size.points.3d, 
                  color = Z.class[, 4])
            }
        }
        for (i in 2:p) graph.coord <- rbind(graph.coord, z.axes[[i]][c(1, 
            nrow(z.axes[[i]])), 1:3])
        open3d()
        text3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
        points3d(graph.coord, size = 0)
        aspect3d(aspect.3d)
        if (!is.null(specify.classes)) 
            .samples.plot.3d(Z = Z, Z.means.mat = Z.means.mat, 
                class.vec = class.vec, specify.classes = specify.classes, 
                size.points.3d = size.points.3d)
        if (ID.labs) 
            text3d(Z[, 1:3], text = ID.3d, col = "black", font = font.3d, 
                cex = cex.3d)
        for (i in 1:p) {
            axis <- z.axes[[i]]
            text3d(axis[axis[, 4] == max(axis[, 4]), 1:3], text = z.axes.names[i], 
                col = col.text.3d, font = font.3d, cex = cex.3d)
            axis <- axis[axis[, 4] != max(axis[, 4]), 1:3]
            lines3d(axis[c(1, nrow(axis)), 1:3], size = size.ax.3d, 
                color = ax.col[i])
        }
        if (plot.class.means) 
            points3d(x = Z.means.mat[, 1], y = Z.means.mat[, 
                2], z = Z.means.mat[, 3], size = size.means.3d, 
                color = Z.means.mat[, 5])
        if (!(is.null(X.new.profiles))) {
            pch.new <- paste("N", 1:nrow(X.new.profiles), sep = "")
            Z.new <- X.new.profiles %*% GH.mat[, 1:3]
            text3d(Z.new, text = pch.new, adj = c(0.5, 0.5), 
                col = "black", font = font.3d, cex = cex.3d)
        }
        rgl.surface(x = c(min(graph.coord[, 1]), max(graph.coord[, 
            1])) * factor.x, z = c(min(graph.coord[, 2]), max(graph.coord[, 
            2])) * factor.y, y = matrix(0, nrow = 2, ncol = 2), 
            coords = c(1, 3, 2), col = col.plane.3d, alpha = alpha.3d)
        axes3d(col = ax.col.3d)
        title3d(main = Titles.3d[1], sub = Titles.3d[2], xlab = Titles.3d[3], 
            ylab = Titles.3d[4], zlab = Titles.3d[5])
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
