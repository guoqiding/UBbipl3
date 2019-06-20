Colour.Pixel <-
function (x, y, mat, colours = 1:nrow(mat), dist = "Pythagoras") 
{
    dist.vec <- apply(mat, 1, function(row, x, y, dist) {
        vec <- c(x, y, rep(0, ncol(mat) - 2))
        if (dist == "Pythagoras") 
            sqrt(sum((vec - row)^2))
    }, x = x, y = y, dist = dist)
    colours[dist.vec == min(dist.vec)]
}
