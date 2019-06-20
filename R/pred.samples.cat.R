pred.samples.cat <-
function (samples.mat, Z.rj, var.j, cij2) 
{
    n <- ncol(Z.rj)
    predictions <- rep(NA, nrow(samples.mat))
    for (i in 1:nrow(samples.mat)) {
        temp.mat <- samples.mat[i, ] %*% matrix(1, nrow = 1, 
            ncol = n) - Z.rj
        dist.vec <- apply(temp.mat, 2, function(a) sum(a * a))
        dist.vec <- dist.vec + cij2
        temp.vec <- dist.vec == min(dist.vec)
        cats <- levels(var.j)
        var.char <- as.character(var.j)
        predicted <- as.vector(unique(cats[match((var.char[temp.vec]), 
            cats)]))
        if (length(predicted) > 1) 
            warning(paste("Point predicted by all of: ", predicted))
        predictions[i] <- predicted[1]
    }
    return(predictions)
}
