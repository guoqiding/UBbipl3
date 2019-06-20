pred.regions.cat <-
function (Z.rj, cij2, var.j, x.grid = 0.05, y.grid = 0.05, plot.symbol = 15, 
    size = 0.25, colours = 1:10) 
{
    usr <- par("usr")
    grid <- expand.grid(seq(from = usr[1], to = usr[2], by = x.grid), 
        seq(from = usr[3], to = usr[4], by = y.grid))
    apply(grid, 1, function(row, plot.symbol, size, Z.rj, cij2, 
        var.j, colours) {
        points(x = row[1], y = row[2], pch = plot.symbol, cex = size, 
            col = Colour.Pixel.cat(x = row[1], y = row[2], Z.rj, 
                var.j, cij2, colours))
    }, plot.symbol = plot.symbol, size = size, Z.rj = Z.rj, cij2 = cij2, 
        var.j = var.j, colours = colours)
}
