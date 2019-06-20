biad.ss <-
function (X) 
{
    row.means <- apply(X, 1, mean)
    col.means <- apply(X, 2, mean)
    p <- nrow(X)
    q <- ncol(X)
    Interaction.mat <- X - matrix(1, nrow = p, ncol = 1) %*% 
        (matrix(1, nrow = 1, ncol = p) %*% X/p) - X %*% (rep(1, 
        q)/q) %*% matrix(1, nrow = 1, ncol = q) + mean(X)
    uit.svd <- svd(Interaction.mat)
    C.mat <- as.matrix(uit.svd$u %*% diag(sqrt(uit.svd$d)))
    dimnames(C.mat) <- list(dimnames(X)[[1]], NULL)
    D.mat <- uit.svd$v %*% diag(sqrt(uit.svd$d))
    dimnames(D.mat) <- list(dimnames(X)[[2]], NULL)
    e.vals <- uit.svd$d
    X.hat <- matrix(1, nrow = p, ncol = 1) %*% (matrix(1, nrow = 1, 
        ncol = p) %*% X/p) + X %*% (rep(1, q)/q) %*% matrix(1, 
        nrow = 1, ncol = q)
    -mean(X) + C.mat[, 1] %*% t(D.mat[, 1]) + C.mat[, 2] %*% 
        t(D.mat[, 2])
    dimnames(X.hat) <- dimnames(X)
    if (p < q) 
        aantal <- p - 1
    else aantal <- q - 1
    Tot.SS <- sum(X * X)
    SS.Mean <- ((sum(X))^2)/(p * q)
    SS.Tot.Meancorrected <- Tot.SS - SS.Mean
    SS.Rows <- sum(apply(X, 1, sum)^2)/q - SS.Mean
    SS.Cols <- sum(apply(X, 2, sum)^2)/p - SS.Mean
    SS.Interaction.Total <- sum(diag(Interaction.mat %*% t(Interaction.mat)))
    SS.Interaction.components <- (diag(t(C.mat) %*% C.mat) * 
        diag(t(D.mat) %*% D.mat))
    Vec.SS <- c(SS.Rows, SS.Cols, SS.Interaction.Total, SS.Mean, 
        Tot.SS)
    Vec.DF <- c(p - 1, q - 1, (p - 1) * (q - 1), 1, p * q)
    Vec.MSS <- Vec.SS/Vec.DF
    Source <- c("Row main effect", "Column main effect", "Interaction", 
        "Mean", "Total (Uncor)")
    out <- data.frame(Source = Source, Sum_of_squares = Vec.SS, 
        DF = Vec.DF, Mean_ss = Vec.MSS)
    SS.Interaction.components <- c(SS.Interaction.Total, SS.Interaction.components[1:aantal])
    SS.Interaction.DF <- c((p - 1) * (q - 1), seq(from = p + 
        q - 1 - 2 * 1, to = p + q - 1 - 2 * aantal, by = -2))
    MSS.Interaction.components <- SS.Interaction.components/SS.Interaction.DF
    Heading <- c("Row:Column Interaction", paste("Multiplicative term ", 
        1:aantal))
    out.2 <- data.frame(Heading, SS.Interaction.components, SS.Interaction.DF, 
        MSS.Interaction.components)
    names(out.2) <- names(out)
    list(X = X, X.hat = X.hat, Interaction.mat = Interaction.mat, 
        C.mat = C.mat, D.mat = D.mat, SS.Tot.Meancorrected = SS.Tot.Meancorrected, 
        SS.Interaction.Total = SS.Interaction.Total, SS.Interaction.components = SS.Interaction.components, 
        cumsum = cumsum(SS.Interaction.components), ANOVA.Table = rbind(out, 
            out.2))
}
