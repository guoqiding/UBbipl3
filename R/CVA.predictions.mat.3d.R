CVA.predictions.mat.3d <-
function (X, B) 
{
    X <- as.matrix(X)
    p <- ncol(X)
    X.centred <- scale(X, center = TRUE, scale = FALSE)
    predictions.for.centredX <- vector("list", p)
    predictions.for.originalX <- vector("list", p)
    reconstr.error <- rep(NA, p)
    for (r in 1:p) {
        Br <- B[, 1:r]
        Brr <- solve(B)[1:r, ]
        predict.centredX <- X.centred %*% Br %*% Brr
        predictions.for.centredX[[r]] <- predict.centredX
        predictions.for.originalX[[r]] <- scale(predict.centredX, 
            center = -attr(X.centred, "scaled:center"), scale = FALSE)
        dimnames(predictions.for.centredX[[r]])[[2]] <- dimnames(X)[[2]]
        dimnames(predictions.for.originalX[[r]])[[2]] <- dimnames(X)[[2]]
        reconstr.error[r] <- sum((X.centred - predict.centredX)^2)
    }
    list(predictions.for.centredX = predictions.for.centredX, 
        predictions.for.originalX = predictions.for.originalX, 
        reconstr.error = reconstr.error)
}
