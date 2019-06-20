projections <-
function (from, to, width, minvec) 
{
    bigmat <- scale(cbind(from, to), center = rep(minvec, 2), 
        scale = FALSE)
    apply(bigmat, 1, function(x, width) {
        my.arrow(x[c(1, 4)], x[c(3, 6)], x[c(2, 5)], width = width, 
            colour = "black", size = 2)
    }, width = width)
}
