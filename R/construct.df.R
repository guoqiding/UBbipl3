construct.df <-
function (dat.list = RSACrime.data, year = 1) 
{
    dat.mat <- dat.list[[1]][, year, drop = FALSE]
    for (i in 2:length(dat.list)) dat.mat <- cbind(dat.mat, dat.list[[i]][, 
        year, drop = FALSE])
    dimnames(dat.mat)[[2]] <- names(dat.list)
    dat.mat
}
