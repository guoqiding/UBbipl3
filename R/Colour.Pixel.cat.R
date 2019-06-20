Colour.Pixel.cat <-
function (x, y, Z.rj, var.j, cij2, colours = 1:length(levels(var.j))) 
{
    point <- matrix(c(x, y), nrow = 2, ncol = 1)
    n <- ncol(Z.rj)
    temp.mat <- point %*% matrix(1, nrow = 1, ncol = n) - Z.rj
    dist.vec <- apply(temp.mat, 2, function(a) sum(a * a))
    dist.vec <- dist.vec + cij2
    temp.vec <- dist.vec == min(dist.vec)
    cats <- levels(var.j)
    var.char <- as.character(var.j)
    predicted <- as.vector(unique(cats[match((var.char[temp.vec]), 
        cats)]))
    if (length(predicted) > 1) 
        warning(paste("Point predicted by all of: ", predicted))
    colours[match((var.char[temp.vec])[1], cats)]
}
