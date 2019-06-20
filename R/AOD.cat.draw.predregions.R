AOD.cat.draw.predregions <-
function (X = Ocotea.cat[, 3:8], group.vec = Ocotea.cat[, 2], 
    exp.factor = 1.2, parplotmar = c(3, 3, 3, 3), dist.metric = c("EMC", 
        "ChiSq"), scaffolding = c(1, 2), plot.CLPs = NULL, predregions = 1, 
    x.grid = 0.015, y.grid = 0.015, pred.region.symbol = 15, 
    pred.region.symbol.size = 0.2, pch.sample = c(22, 24, 25), 
    pch.group = c(22, 24, 25), pch.sample.lwd = 2, legend.plot = TRUE, 
    CLP.LabelsInPredRegs = FALSE, region.colours = UBcolours[1:5], 
    sample.labels = TRUE, sample.labels.cex = 0.7, text.size = c(1, 
        1), colours = UBcolours, cij.null = FALSE, predict.samples = NULL, 
    labels = NULL, CLP.multiplier = 1, jitter.factor = 0) 
{
    if (length(scaffolding) != 2) 
        stop("Scaffolding must be a two-componnet integer-valued vector \n")
    predictions.made <- NULL
    Pythagoras.dist <- function(Xmat) {
        n <- nrow(Xmat)
        B <- Xmat %*% t(Xmat)
        Delta <- matrix(0, nrow = n, ncol = n)
        for (i in 1:n) for (j in 1:n) Delta[i, j] <- B[i, i] + 
            B[j, j] - 2 * B[i, j]
        Delta
    }
    dist.metric <- dist.metric[1]
    n <- nrow(X)
    p <- ncol(X)
    Ggroupmat <- indmat(group.vec)
    Nmat <- t(Ggroupmat) %*% Ggroupmat
    nvec <- diag(Nmat)
    K <- length(nvec)
    Glist <- vector("list", p)
    for (i in 1:p) Glist[[i]] <- indmat(X[, i])
    Gmat.colnames <- NULL
    for (i in 1:p) Gmat.colnames <- append(Gmat.colnames, dimnames(Glist[[i]])[[2]])
    Gmat <- matrix(unlist(Glist), nrow = n)
    ncat <- ncol(Gmat)
    dimnames(Gmat) <- list(dimnames(X)[[1]], Gmat.colnames)
    Llist <- vector("list", p)
    for (j in 1:p) Llist[[j]] <- as.matrix(t(Glist[[j]]) %*% 
        Glist[[j]])
    L.mat <- diag(diag(t(Gmat) %*% Gmat))
    Dmat.list <- vector("list", p)
    if (dist.metric == "EMC") {
        for (j in 1:p) Dmat.list[[j]] <- -(0.5) * (matrix(1, 
            nrow = n, ncol = n) - Glist[[j]] %*% t(Glist[[j]]))
        Dmat <- -(0.5) * (p * matrix(1, nrow = n, ncol = n) - 
            Gmat %*% t(Gmat))
    }
    if (dist.metric == "ChiSq") {
        for (j in 1:p) Dmat.list[[j]] <- -(0.5) * (Pythagoras.dist((p^(-0.5) * 
            Glist[[j]] %*% solve(sqrt(Llist[[j]])))))
        Dmat <- -(0.5) * (Pythagoras.dist((p^(-0.5) * Gmat %*% 
            solve(sqrt(L.mat)))))
    }
    Sum.Djs <- Dmat.list[[1]]
    for (j in 2:p) Sum.Djs <- Sum.Djs + Dmat.list[[j]]
    Centmat.n <- diag(n) - matrix(1, nrow = n, ncol = n)/n
    B.mat <- Centmat.n %*% Dmat %*% Centmat.n
    Bmat.list <- vector("list", p)
    for (j in 1:p) Bmat.list[[j]] <- Centmat.n %*% Dmat.list[[j]] %*% 
        Centmat.n
    Sum.Bjs <- Bmat.list[[1]]
    for (j in 2:p) Sum.Bjs <- Sum.Bjs + Bmat.list[[j]]
    PCO.B <- svd(B.mat)
    lambda.mat <- zapsmall(diag(PCO.B$d))
    lambda.mat.inverse <- ifelse(lambda.mat > 0, 1/lambda.mat, 
        0)
    Y.mat <- PCO.B$v %*% sqrt(lambda.mat)
    Zmat.list <- vector("list", p)
    for (j in 1:p) Zmat.list[[j]] <- lambda.mat.inverse %*% t(Y.mat) %*% 
        Bmat.list[[j]]
    Dbar <- solve(Nmat) %*% t(Ggroupmat) %*% Dmat %*% Ggroupmat %*% 
        solve(Nmat)
    if (nrow(Dbar) == 1) 
        Deltamat <- Dbar - 0.5 * (Dbar %*% matrix(1, nrow = K, 
            ncol = K) + matrix(1, nrow = K, ncol = K) %*% Dbar)
    else Deltamat <- Dbar - 0.5 * (diag(diag(Dbar)) %*% matrix(1, 
        nrow = K, ncol = K) + matrix(1, nrow = K, ncol = K) %*% 
        diag(diag(Dbar)))
    Centmat.K <- diag(K) - matrix(1, nrow = K, ncol = K)/K
    T.ss <- -sum(Dmat)/sum(nvec)
    B.ss <- -(t(nvec) %*% Deltamat %*% nvec)/sum(nvec)
    W.ss <- T.ss - B.ss
    B.delta <- Centmat.K %*% Deltamat %*% Centmat.K
    PCO.B.delta <- svd(B.delta)
    lambda.delta <- zapsmall(diag(PCO.B.delta$d))
    lambda.delta.inverse <- ifelse(lambda.delta > 0, 1/lambda.delta, 
        0)
    Y.delta <- PCO.B.delta$v %*% sqrt(lambda.delta)
    Y.delta <- PCO.B.delta$u %*% sqrt(lambda.delta)
    par(pty = "s", mar = parplotmar)
    plot(Y.delta[, scaffolding], asp = 1)
    points(x = 0, y = 0, pch = "O", col = "red", cex = 2)
    translation.vector <- matrix(nvec/n, nrow = K, ncol = 1)
    Y.delta.transl <- t(translation.vector) %*% Y.delta
    points(x = Y.delta.transl[1], y = Y.delta.transl[2], pch = "G", 
        col = "blue", cex = 2)
    Y.delta.G <- t(Y.delta) - t(Y.delta) %*% translation.vector %*% 
        matrix(1, nrow = 1, ncol = K)
    windows()
    par(pty = "s", mar = parplotmar)
    plot(t(Y.delta.G)[, scaffolding], asp = 1)
    points(x = 0, y = 0, pch = "G", col = "red", cex = 2)
    Y.samples <- lambda.delta.inverse %*% t(Y.delta) %*% sweep(t(Ggroupmat), 
        1, nvec, "/") %*% B.mat
    windows()
    par(pty = "s", mar = parplotmar)
    plot(t(Y.samples)[, scaffolding], asp = 1, pch = 15, cex = 0.8)
    points(t(Y.delta.G)[, scaffolding], asp = 1, pch = 16, col = "red", 
        cex = 1.25)
    if (sample.labels) 
        text(t(Y.samples)[, scaffolding], labels = rownames(X), 
            cex = 0.8, pos = 3)
    CLP.list <- vector("list", p)
    for (j in 1:p) {
        CLP.list[[j]] <- lambda.delta.inverse %*% t(Y.delta) %*% 
            sweep(t(Ggroupmat), 1, nvec, "/") %*% Bmat.list[[j]]
        colnames(CLP.list[[j]]) <- X[, j]
    }
    Zmat <- t(Y.samples)[, scaffolding]
    windows()
    par(pty = "s", mar = parplotmar)
    limx <- range(Zmat[, 1])
    limy <- range(Zmat[, 2])
    plot(x = Zmat[, 1], y = Zmat[, 2], xlim = limx * exp.factor, 
        ylim = limy * exp.factor, type = "n", asp = 1, xlab = "", 
        ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
        main = "AOD biplot: categorical vars")
    points(t(Y.delta.G)[, scaffolding], pch = pch.group, col = UBcolours[1:K], 
        cex = 1.5)
    points(x = 0, y = 0, pch = "G", col = "black", cex = 1.2)
    for (k in 1:K) {
        select.grp <- as.logical(Ggroupmat[, k])
        points(x = jitter(Zmat[select.grp, 1], factor = 5), y = jitter(Zmat[select.grp, 
            2], factor = 5), col = UBcolours[k], pch = pch.sample[k])
        if (sample.labels) 
            text(x = Zmat[select.grp, 1], y = Zmat[select.grp, 
                2], labels = rownames(X)[select.grp], cex = 0.8, 
                pos = 3)
    }
    names.grps <- colnames(Ggroupmat)
    if (legend.plot) 
        legend(x = "bottomright", legend = names.grps, pch = pch.group, 
            col = colours[1:K], cex = text.size[[1]], bg = "white")
    if (!is.null(plot.CLPs)) 
        for (i in plot.CLPs) {
            select.CLP <- t(CLP.list[[i]])[, scaffolding]
            CLPlabel <- NULL
            for (sampl in 1:nrow(select.CLP)) {
                name <- as.character(X[sampl, i])
                if (!(name %in% CLPlabel)) {
                  CLPlabel <- append(CLPlabel, name)
                  points(x = select.CLP[sampl, 1] * CLP.multiplier, 
                    y = select.CLP[sampl, 2] * CLP.multiplier, 
                    pch = 24, col = colours[i], cex = 1)
                  text(x = select.CLP[sampl, 1] * CLP.multiplier, 
                    y = select.CLP[sampl, 2] * CLP.multiplier, 
                    labels = CLPlabel[length(CLPlabel)], col = colours[i], 
                    pos = 3, cex = 0.8)
                }
            }
        }
    if (!is.null(predregions)) {
        windows()
        for (j in predregions) {
            par(pty = "s", mar = parplotmar)
            plot(x = Zmat[, 1], y = Zmat[, 2], xlim = limx * 
                exp.factor, ylim = limy * exp.factor, type = "n", 
                asp = 1, xlab = "", ylab = "", xaxs = "i", xaxt = "n", 
                yaxt = "n", yaxs = "i", main = paste("Prediction regions for variable ", 
                  j, sep = "Var "))
            var.j <- X[, j]
            select.var <- unique(round(t(CLP.list[[j]]), digits = 8))
            Z.rj <- CLP.list[[j]][scaffolding, ]
            D.j <- Dmat.list[[j]]
            bar.D.j <- sum(D.j)/(n * n)
            dist2.P0fromG <- bar.D.j * rep(1, n) - 2 * matrix(1, 
                nrow = 1, ncol = n) %*% D.j %*% diag(n)/n
            if (!cij.null) 
                cij2 <- dist2.P0fromG - apply(Z.rj, 2, function(x) sum(x * 
                  x))
            else cij2 <- 0
            pred.regions.cat(Z.rj = Z.rj, cij2 = cij2, x.grid = x.grid, 
                y.grid = y.grid, var.j = var.j, plot.symbol = pred.region.symbol, 
                size = pred.region.symbol.size, colours = region.colours)
            number <- readline("Enter Number of regions in biplot \n")
            sides <- readline("Enter Number of corners and then select the corner points \n")
            pol.points <- locator(sides)
            text(pol.points, label = 1:sides, xpd = NA)
            windows()
            par(pty = "s", mar = parplotmar)
            plot(x = Zmat[, 1], y = Zmat[, 2], xlim = limx * 
                exp.factor, ylim = limy * exp.factor, type = "n", 
                asp = 1, xlab = "", ylab = "", xaxs = "i", xaxt = "n", 
                yaxt = "n", yaxs = "i", main = paste("Prediction regions for variable Var", 
                  j, "cij.null =", cij.null, sep = " "))
            for (ask in 1:number) {
                plot.corners <- readline("Specify corners of polygon plot as a vector \n")
                plot.corners <- eval(parse(text = plot.corners))
                region.col <- readline("Specify colour of polygon plot without using quotes  \n")
                polygon(x = pol.points$x[plot.corners], y = pol.points$y[plot.corners], 
                  col = region.col, border = region.col)
            }
            points(t(Y.delta.G)[, scaffolding], pch = pch.group, 
                col = "black", bg = "black", cex = 1.25)
            text(t(Y.delta.G)[, scaffolding], labels = names.grps, 
                cex = sample.labels.cex, pos = 3)
            samples.cols <- region.colours[unclass(X[, predregions])]
            for (k in 1:K) {
                select.group <- as.logical(Ggroupmat[, k])
                points(x = jitter(t(Y.samples)[select.group, 
                  1], factor = 5), y = jitter(t(Y.samples)[select.group, 
                  2], factor = 5), col = "black", bg = samples.cols[select.group], 
                  pch = pch.sample[k], lwd = 2, cex = pch.sample.lwd)
                text(x = t(Y.samples)[select.group, 1], y = t(Y.samples)[select.group, 
                  2], labels = rownames(X)[select.group], cex = sample.labels.cex, 
                  pos = 3)
            }
            points(0, 0, pch = 3)
            if (CLP.LabelsInPredRegs) {
                points(x = select.var[, 1], y = select.var[, 
                  2], pch = 16, cex = 0.6, col = "black")
                text(x = select.var[, 1], y = select.var[, 2], 
                  label = rownames(select.var), pos = 1, offset = 0.1, 
                  cex = 0.8)
            }
            names <- levels(var.j)
            legend(x = "bottomleft", legend = names, pch = rep(15, 
                length(names)), col = region.colours, cex = text.size[[1]], 
                bg = "white")
            legend(x = "bottomright", legend = names.grps, pch = pch.group, 
                col = c("black", bg = "black"), cex = text.size[[1]], 
                bg = "white")
            windows()
        }
    }
    if (!is.null(predict.samples)) {
        samples.mat <- t(Y.samples)[predict.samples, scaffolding, 
            drop = FALSE]
        predictions.made <- rep(NA, nrow(samples.mat))
        for (j in 1:p) {
            var.j <- X[, j]
            select.var <- t(CLP.list[[j]])
            Z.rj <- CLP.list[[j]][scaffolding, ]
            D.j <- Dmat.list[[j]]
            bar.D.j <- sum(D.j)/(n * n)
            dist2.P0fromG <- bar.D.j * rep(1, n) - 2 * matrix(1, 
                nrow = 1, ncol = n) %*% D.j %*% diag(n)/n
            if (!cij.null) 
                cij2 <- dist2.P0fromG - apply(Z.rj, 2, function(x) sum(x * 
                  x))
            else cij2 <- 0
            pred <- pred.samples.cat(samples.mat = samples.mat, 
                Z.rj = Z.rj, var.j = var.j, cij2 = cij2)
            predictions.made <- data.frame(predictions.made, 
                pred)
        }
        predictions.made <- predictions.made[, -1]
        colnames(predictions.made) <- colnames(X)
        rownames(predictions.made) <- rownames(X)[predict.samples]
    }
    return(predictions.made)
}
