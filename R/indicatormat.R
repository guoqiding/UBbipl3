indicatormat <-
function (table) 
{
    nr <- nrow(table)
    nc <- ncol(table)
    fvec <- c(t(table))
    tmp1 <- diag(nr)
    tmp2 <- diag(nc)
    tmpmat <- NULL
    for (i in 1:nr) for (j in 1:nc) tmpmat <- rbind(tmpmat, cbind(tmp1[i, 
        , drop = F], tmp2[j, , drop = F]))
    Indicatormat <- tmpmat[rep((1:nrow(tmpmat)), fvec), ]
    dimnames(Indicatormat) <- list(NULL, c(dimnames(table)[[1]], 
        dimnames(table)[[2]]))
    Indicatormat
}
