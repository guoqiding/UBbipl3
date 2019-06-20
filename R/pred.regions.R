pred.regions <-
function (mat, x.grid = 0.05, y.grid = 0.05, plot.symbol = 15, 
    size = 0.25, colours = 1:10) 
{
    usr <- par("usr")
    grid <- expand.grid(seq(from = usr[1], to = usr[2], by = x.grid), 
        seq(from = usr[3], to = usr[4], by = y.grid))
    apply(grid, 1, function(row, plot.symbol, size, mat, colours) points(x = row[1], 
        y = row[2], pch = plot.symbol, cex = size, col = Colour.Pixel(x = row[1], 
            y = row[2], mat, colours)), plot.symbol = plot.symbol, 
        size = size, mat = mat, colours = colours)
}
