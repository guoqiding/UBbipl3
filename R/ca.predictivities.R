ca.predictivities <-
function (data = as.matrix(SACrime08.data), out = 1:6) 
{
    nc <- ncol(data)
    nr <- nrow(data)
    nm <- min(nc, nr)
    R.mat <- diag(apply(data, 1, sum))
    C.mat <- diag(apply(data, 2, sum))
    E.mat <- R.mat %*% matrix(1, nrow = nr, ncol = nc) %*% C.mat/sum(data)
    RHalf <- sqrt(R.mat)
    CHalf <- sqrt(C.mat)
    RMinHalf <- solve(RHalf)
    CMinHalf <- solve(CHalf)
    dev.mat <- data - E.mat
    weighted.dev.mat <- (sqrt(solve(R.mat)) %*% dev.mat %*% sqrt(solve(C.mat)))
    svd.weighted.dev.mat <- svd(weighted.dev.mat)
    sing.values <- svd.weighted.dev.mat$d
    Sigma2 <- diag(sing.values)^2
    Quality <- (100 * cumsum(diag(Sigma2)))/sum(diag(Sigma2))
    names(Quality) <- paste("Dim", 1:length(Quality))
    NoemerMat <- diag(diag(svd.weighted.dev.mat$v %*% Sigma2 %*% 
        t(svd.weighted.dev.mat$v)))
    Weights <- diag(NoemerMat)/sum(Sigma2)
    Uit.list.cols <- vector("list", nc)
    Uit.list.Adeq <- vector("list", nc)
    Uit.list.Xhat <- vector("list", nm)
    for (i in 1:nm) {
        Jr <- matrix(0, nrow = nm, ncol = nm)
        Jr[1:i, 1:i] <- diag(i)
        Uit.list.cols[[i]] <- diag(diag(svd.weighted.dev.mat$v %*% 
            Sigma2 %*% Jr %*% t(svd.weighted.dev.mat$v))) %*% 
            solve(NoemerMat)
        Uit.list.Adeq[[i]] <- diag(diag(svd.weighted.dev.mat$v %*% 
            Jr %*% t(svd.weighted.dev.mat$v)))
        Uit.list.Xhat[[i]] <- svd.weighted.dev.mat$u %*% diag(svd.weighted.dev.mat$d) %*% 
            Jr %*% t(svd.weighted.dev.mat$v)
        dimnames(Uit.list.Xhat[[i]]) <- dimnames(data)
    }
    Adequacies <- round(diag(Uit.list.Adeq[[1]]), digits = 4)
    for (i in 2:nm) Adequacies <- rbind(Adequacies, round(diag(Uit.list.Adeq[[i]]), 
        digits = 4))
    dimnames(Adequacies) <- list(paste("Dim", 1:nm, sep = "_"), 
        dimnames(data)[[2]])
    Uit.list.samples <- vector("list", nc)
    for (i in 1:nm) {
        Jr <- matrix(0, nrow = nm, ncol = nm)
        Jr[1:i, 1:i] <- diag(i)
        Uit.list.samples[[i]] <- diag(diag(svd.weighted.dev.mat$u %*% 
            Sigma2 %*% Jr %*% t(svd.weighted.dev.mat$u))) %*% 
            solve(diag(diag(svd.weighted.dev.mat$u %*% Sigma2 %*% 
                t(svd.weighted.dev.mat$u))))
    }
    Axis.predictivities <- round(diag(Uit.list.cols[[1]]), digits = 4)
    for (i in 2:nm) Axis.predictivities <- rbind(Axis.predictivities, 
        round(diag(Uit.list.cols[[i]]), digits = 4))
    dimnames(Axis.predictivities) <- list(paste("Dim", 1:nm, 
        sep = "_"), dimnames(data)[[2]])
    Sample.predictivities <- round(diag(Uit.list.samples[[1]]), 
        digits = 4)
    for (i in 2:nm) Sample.predictivities <- rbind(Sample.predictivities, 
        round(diag(Uit.list.samples[[i]]), digits = 4))
    dimnames(Sample.predictivities) <- list(paste("Dim", 1:nm, 
        sep = "_"), dimnames(data)[[1]])
    out.list <- list(Quality = round(Quality, digits = 4), Weights = Weights, 
        Adequacies = Adequacies, Axis.predictivities = Axis.predictivities, 
        Sample.predictivities = Sample.predictivities, X.hat = Uit.list.Xhat)
    out.list[out]
}
