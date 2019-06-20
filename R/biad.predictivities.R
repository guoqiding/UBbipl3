biad.predictivities <-
function (X = wheat.data, e.vects = 1:ncol(X), X.new.rows = NULL, 
    X.new.columns = NULL, add.maineffect = FALSE, biad.variant = c("Xmat", 
        "XminMeanMat", "InteractionMat"), predictions.dim = c("All", 
        1:min(nrow(X), ncol(X)))) 
{
    if (!is.null(X.new.rows)) 
        if (!is.matrix(X.new.rows) | !identical(ncol(X), ncol(X.new.rows))) 
            stop("X.new.rows must be a matrix of size k x q \n")
    if (!is.null(X.new.columns)) 
        if (!is.matrix(X.new.columns) | !identical(nrow(X), nrow(X.new.columns))) 
            stop("X.new.columns must be a matrix of size p x m \n")
    New.Row.predictivities <- NULL
    New.Column.predictivities <- NULL
    biad.variant <- biad.variant[1]
    predictions.dim <- predictions.dim[1]
    Xmean <- mean(X)
    p <- nrow(X)
    q <- ncol(X)
    r <- min(p, q)
    Interaction.mat <- X - matrix(1, nrow = p, ncol = 1) %*% 
        (matrix(1, nrow = 1, ncol = p) %*% X/p) - X %*% (rep(1, 
        q)/q) %*% matrix(1, nrow = 1, ncol = q) + mean(X)
    XMinOverallMean <- X - Xmean
    main.effects.rows <- apply(X, 1, mean) - mean(X)
    main.effects.cols <- apply(X, 2, mean) - mean(X)
    Residual.mat <- X - matrix(apply(X, 1, mean), ncol = 1) %*% 
        matrix(1, nrow = 1, ncol = ncol(X)) - matrix(1, nrow = nrow(X), 
        ncol = 1) %*% matrix(apply(X, 2, mean), nrow = 1) + mean(X)
    svd.X <- svd(X)
    svd.XMinMean <- svd(XMinOverallMean)
    svd.InteractMat <- svd(Residual.mat)
    U.full.X <- svd(X %*% t(X))$u
    U.full.XMinMean <- svd(XMinOverallMean %*% t(XMinOverallMean))$u
    U.full.InteractMat <- svd(Residual.mat %*% t(Residual.mat))$u
    if (biad.variant == "Xmat") {
        Umat <- svd.X$u
        Sigma <- diag(svd.X$d)
        Vmat <- svd.X$v
        SVDmat <- X
        U.full <- U.full.X
        if (!is.null(X.new.rows)) 
            H.new.rows <- X.new.rows
        if (!is.null(X.new.columns)) 
            H.new.columns <- X.new.columns
    }
    if (biad.variant == "XminMeanMat") {
        Umat <- svd.XMinMean$u
        Sigma <- diag(svd.XMinMean$d)
        Vmat <- svd.XMinMean$v
        SVDmat <- XMinOverallMean
        U.full <- U.full.XMinMean
        if (!is.null(X.new.rows)) 
            H.new.rows <- X.new.rows - matrix(1, nrow = nrow(X.new.rows), 
                ncol = ncol(X.new.rows)) * Xmean
        if (!is.null(X.new.columns)) 
            H.new.columns <- X.new.columns - matrix(1, nrow = nrow(X.new.columns), 
                ncol = ncol(X.new.columns)) * Xmean
    }
    if (biad.variant == "InteractionMat") {
        Umat <- svd.InteractMat$u
        Sigma <- diag(svd.InteractMat$d)
        Vmat <- svd.InteractMat$v
        SVDmat <- Residual.mat
        U.full <- U.full.InteractMat
        if (!is.null(X.new.rows)) {
            X.new.rows.row.means <- apply(X.new.rows, 1, mean)
            X.column.means <- apply(X, 2, mean)
            H.new.rows <- X.new.rows - matrix(X.new.rows.row.means, 
                byrow = FALSE, nrow = nrow(X.new.rows), ncol = ncol(X)) - 
                matrix(X.column.means, byrow = TRUE, nrow = nrow(X.new.rows), 
                  ncol = ncol(X)) + matrix(1, nrow = nrow(X.new.rows), 
                ncol = ncol(X)) * Xmean
        }
        if (!is.null(X.new.columns)) {
            X.new.columns.column.means <- apply(X.new.columns, 
                2, mean)
            X.row.means <- apply(X, 1, mean)
            H.new.columns <- X.new.columns - matrix(X.row.means, 
                byrow = FALSE, nrow = nrow(X), ncol = ncol(X.new.columns)) - 
                matrix(X.new.columns.column.means, byrow = TRUE, 
                  nrow = nrow(X), ncol = ncol(X.new.columns)) + 
                matrix(1, nrow = nrow(X), ncol = ncol(X.new.columns)) * 
                  Xmean
        }
    }
    sing.values <- diag(Sigma)
    Quality <- cumsum((sing.values^2)[e.vects][1:r])/sum((sing.values^2)) * 
        100
    Sigma2 <- Sigma %*% Sigma
    NumMat <- diag(diag(Vmat %*% Sigma2 %*% t(Vmat)))
    Uit.list.cols <- vector("list", q)
    Uit.list.col.adequacies <- vector("list", q)
    Uit.list.row.adequacies <- vector("list", p)
    Uit.list.Xhat <- vector("list", min(p, q))
    for (i in 1:min(p, q)) {
        Jr <- matrix(0, nrow = min(p, q), ncol = min(p, q))
        Jr[1:i, 1:i] <- diag(i)
        Uit.list.cols[[i]] <- diag(diag(Vmat %*% Sigma2 %*% Jr %*% 
            t(Vmat))) %*% solve(NumMat)
        Uit.list.col.adequacies[[i]] <- diag(diag(Vmat %*% Jr %*% 
            t(Vmat)))
        Uit.list.Xhat[[i]] <- Umat %*% Sigma %*% Jr %*% t(Vmat)
        if (biad.variant == "InteractionMat") {
            if (add.maineffect) 
                Uit.list.Xhat[[i]] <- scale(Uit.list.Xhat[[i]], 
                  scale = FALSE, center = -main.effects.cols)
        }
        dimnames(Uit.list.Xhat[[i]]) <- dimnames(X)
    }
    for (i in 1:p) {
        Jr <- matrix(0, nrow = p, ncol = p)
        Jr[1:i, 1:i] <- diag(i)
        Uit.list.row.adequacies[[i]] <- diag(diag(U.full %*% 
            Jr %*% t(U.full)))
    }
    Row.adequacies <- round(diag(Uit.list.row.adequacies[[1]]), 
        digits = 4)
    for (i in 2:p) Row.adequacies <- rbind(Row.adequacies, round(diag(Uit.list.row.adequacies[[i]]), 
        digits = 4))
    dimnames(Row.adequacies) <- list(paste("Dim", 1:p, sep = "_"), 
        dimnames(X)[[1]])
    Uit.list.rows <- vector("list", p)
    for (i in 1:min(p, q)) {
        Jr <- matrix(0, nrow = min(p, q), ncol = min(p, q))
        Jr[1:i, 1:i] <- diag(i)
        Xhat <- SVDmat %*% Vmat %*% Jr %*% t(Vmat)
        HtH <- diag(SVDmat %*% t(SVDmat))
        Denom <- diag(1/diag(HtH))
        Uit.list.rows[[i]] <- diag(diag((Xhat) %*% t(Xhat)) * 
            Denom)
    }
    Column.adequacies <- round(diag(Uit.list.col.adequacies[[1]]), 
        digits = 4)
    for (i in 2:min(p, q)) Column.adequacies <- rbind(Column.adequacies, 
        round(diag(Uit.list.col.adequacies[[i]]), digits = 4))
    dimnames(Column.adequacies) <- list(paste("Dim", 1:min(p, 
        q), sep = "_"), dimnames(X)[[2]])
    Column.predictivities <- round(diag(Uit.list.cols[[1]]), 
        digits = 4)
    for (i in 2:min(p, q)) Column.predictivities <- rbind(Column.predictivities, 
        round(diag(Uit.list.cols[[i]]), digits = 4))
    dimnames(Column.predictivities) <- list(paste("Dim", 1:min(p, 
        q), sep = "_"), dimnames(X)[[2]])
    Row.predictivities <- round(diag(Uit.list.rows[[1]]), digits = 4)
    for (i in 2:min(p, q)) Row.predictivities <- rbind(Row.predictivities, 
        round(diag(Uit.list.rows[[i]]), digits = 4))
    dimnames(Row.predictivities) <- list(paste("Dim", 1:min(p, 
        q), sep = "_"), dimnames(X)[[1]])
    if (predictions.dim == "All" | predictions.dim == "all") 
        predictions.rows <- Uit.list.Xhat
    else predictions.rows <- Uit.list.Xhat[[predictions.dim]]
    if (!is.null(X.new.rows)) {
        Uit.list.new.rows <- vector("list", p)
        for (i in 1:min(p, q)) {
            Jr <- matrix(0, nrow = min(p, q), ncol = min(p, q))
            Jr[1:i, 1:i] <- diag(i)
            Numer <- diag(H.new.rows %*% Vmat %*% Jr %*% t(Vmat) %*% 
                t(H.new.rows))
            Denom <- diag(H.new.rows %*% t(H.new.rows))
            Uit.list.new.rows[[i]] <- diag(Numer/Denom)
        }
        New.Row.predictivities <- round(diag(Uit.list.new.rows[[1]]), 
            digits = 4)
        for (i in 2:min(p, q)) New.Row.predictivities <- rbind(New.Row.predictivities, 
            round(diag(Uit.list.new.rows[[i]]), digits = 4))
        dimnames(New.Row.predictivities) <- list(paste("Dim", 
            1:min(p, q), sep = "_"), dimnames(X.new.rows)[[1]])
    }
    if (!is.null(X.new.columns)) {
        Uit.list.new.columns <- vector("list", p)
        for (i in 1:min(p, q)) {
            Jr <- matrix(0, nrow = min(p, q), ncol = min(p, q))
            Jr[1:i, 1:i] <- diag(i)
            Numer <- diag(t(H.new.columns) %*% Umat %*% Jr %*% 
                t(Umat) %*% H.new.columns)
            Denom <- diag(t(H.new.columns) %*% H.new.columns)
            Uit.list.new.columns[[i]] <- diag(Numer/Denom)
        }
        New.Column.predictivities <- round(diag(Uit.list.new.columns[[1]]), 
            digits = 4)
        for (i in 2:min(p, q)) New.Column.predictivities <- rbind(New.Column.predictivities, 
            round(diag(Uit.list.new.columns[[i]]), digits = 4))
        dimnames(New.Column.predictivities) <- list(paste("Dim", 
            1:min(p, q), sep = "_"), dimnames(X.new.columns)[[2]])
    }
    list(Quality = Quality, Column.adequacies = Column.adequacies, 
        Column.predictivities = Column.predictivities, Row.adequacies = Row.adequacies, 
        Row.predictivities = Row.predictivities, main.effects.cols = main.effects.cols, 
        main.effects.rows = main.effects.rows, predictions.rows = predictions.rows, 
        New.Row.predictivities = New.Row.predictivities, New.Column.predictivities = New.Column.predictivities)
}
