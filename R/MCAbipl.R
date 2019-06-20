MCAbipl <-
function (X = MCA.Table.1.data[, -1], e.vects = 1:2, mca.variant = c("indicator", 
    "Burt", "EMC"), category.labels = TRUE, column.points.col = c("red", 
    "red", "blue", "blue", "blue", "blue", "pink", "pink", "pink", 
    "brown", "brown", "brown", "black", "black", "black"), column.points.size = 1.2, 
    column.points.label.size = 0.6, column.points.text = TRUE, 
    constant = 0.05, dim.biplot = c(2, 1, 3), exp.factor = 1.2, 
    main = dimnames(X)[[2]], offset.m = list(rep(-0.1, nrow(X)), 
        rep(-0.1, ncol(X))), ort.lty = 1, parplotmar = rep(1, 
        4), pch.row.points = 16, pch.col.points = 15, pos.m = list(rep(1, 
        nrow(X)), rep(1, ncol(X))), prediction.regions = TRUE, 
    pred.region.symbol = 15, pred.region.symbol.size = 0.3, propshift = 0, 
    reflect = c(FALSE, "x", "y"), region.colours = c("lightgrey", 
        "lightblue", "pink", "cyan"), rotate.degrees = 0, row.points.col = "green", 
    row.points.size = 1.2, row.points.label.size = 0.6, sample.labels = TRUE, 
    text.size = c(0.65, 0.65), Title = "", x.grid = 0.01, y.grid = 0.01, 
    zoomval = NULL, ...) 
{
    mca.variant <- mca.variant[1]
    CLPlist <- NULL
    if (is.na(match(mca.variant, c("indicator", "Burt", "EMC")))) 
        stop("mca variant must be one of: indicator, Burt, EMC \n")
    dim.biplot <- dim.biplot[1]
    if (dim.biplot != 1) {
        constant <- 0
        propshift <- 0
    }
    if (!(dim.biplot == 1 | dim.biplot == 2 | dim.biplot == 3)) 
        stop("Argument dim.biplot must be set to 1 or 2 or 3 \n")
    dim.biplot <- 2
    cat("NOTE we consider a TWO-dimensional biplot \n")
    predictions <- NULL
    p <- ncol(X)
    n <- nrow(X)
    Glist <- vector("list", p)
    for (i in 1:p) Glist[[i]] <- indmat(X[, i])
    Gmat.colnames <- NULL
    for (i in 1:p) Gmat.colnames <- append(Gmat.colnames, dimnames(Glist[[i]])[[2]])
    Gmat <- matrix(unlist(Glist), nrow = n)
    ncat <- ncol(Gmat)
    dimnames(Gmat) <- list(dimnames(X)[[1]], Gmat.colnames)
    Lmat <- diag(apply(Gmat, 2, sum))
    LmatMinHalf <- solve(sqrt(Lmat))
    LmatMinOne <- solve(Lmat)
    BurtMat <- t(Gmat) %*% Gmat
    BurtMat.norm <- LmatMinHalf %*% BurtMat %*% LmatMinHalf/p
    dimnames(BurtMat.norm) <- dimnames(BurtMat)
    EMC.prop.match <- Gmat %*% t(Gmat)/p
    EMC.prop.dissim <- matrix(1, nrow = n, ncol = n) - EMC.prop.match
    svd.pMinhalfGmatLMinhalf <- svd(Gmat %*% LmatMinHalf/sqrt(p), 
        nu = max(n, ncat), nv = max(n, ncat))
    U.mat.full <- svd.pMinhalfGmatLMinhalf$u
    V.mat.full <- svd.pMinhalfGmatLMinhalf$v
    Sigma.mat.full <- diag(svd.pMinhalfGmatLMinhalf$d)
    svd.norm.burt <- svd(BurtMat.norm)
    Vmat <- V.mat.full[, -1]
    Umat <- U.mat.full[, -1]
    SigmaMat <- diag(svd.pMinhalfGmatLMinhalf$d[-1])
    SigmaSquareMat <- diag(svd.norm.burt$d[-1])
    if (dim.biplot == 2) {
        reflect <- reflect[1]
        radns <- pi * rotate.degrees/180
        if (reflect == FALSE) 
            reflect.mat <- diag(2)
        else {
            if (reflect == "x") 
                reflect.mat <- diag(c(1, -1))
            else {
                if (reflect == "y") 
                  reflect.mat <- diag(c(-1, 1))
                else stop("Argument reflect can only set to be NULL, x or y. \n")
            }
        }
    }
    else reflect.mat <- diag(dim.biplot)
    if ((reflect == "x") | (reflect == "y")) 
        if (dim.biplot == 1 && reflect == "y") 
            reflect.mat <- matrix(-1, nrow = 1, ncol = 1)
    if (dim.biplot == 2) 
        rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
            cos(radns)), ncol = 2)
    else rotate.mat <- diag(dim.biplot)
    if (mca.variant == "indicator") {
        par(pty = "s", mar = parplotmar)
        if (ncol(Umat) > ncol(Vmat)) 
            Umat <- Umat[, 1:ncol(Vmat)]
        Z.0 <- Umat %*% SigmaMat[, e.vects]
        Z <- LmatMinHalf %*% Vmat[, e.vects]/sqrt(p)
        Z[, 1:2] <- Z[, 1:2] %*% rotate.mat %*% reflect.mat
        Z.0[, 1:2] <- Z.0[, 1:2] %*% rotate.mat %*% reflect.mat
        plotdata <- rbind(Z[, 1:2], Z.0[, 1:2])
        lims <- c(min(plotdata), max(plotdata))
        plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        points(Z.0, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        points(Z, pch = pch.col.points, col = column.points.col, 
            cex = column.points.size)
        points(0, 0, pch = 3)
        if (sample.labels) 
            text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                offset = offset.m[[1]], cex = text.size[[1]])
        if (category.labels) 
            text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                offset = offset.m[[2]], cex = text.size[[2]])
        if (!is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval1 <- zoom(zoomval)
            plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
                xlim = zoomval1[1:2], ylim = zoomval1[3:4], xlab = "", 
                ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                yaxs = "i", main = Title)
            points(Z.0, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            points(Z, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size)
            points(0, 0, pch = 3)
            if (sample.labels) 
                text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                  offset = offset.m[[1]], cex = text.size[[1]])
            if (category.labels) 
                text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                  offset = offset.m[[2]], cex = text.size[[2]])
        }
        dev.new()
        par(pty = "s", mar = parplotmar)
        Z.0 <- Z.0/p
        plotdata <- rbind(Z[, 1:2], Z.0[, 1:2])
        lims <- c(min(plotdata), max(plotdata))
        plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        points(Z.0, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        points(Z, pch = 15, col = column.points.col, cex = column.points.size)
        points(0, 0, pch = 3)
        if (sample.labels) 
            text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                offset = offset.m[[1]], cex = text.size[[1]])
        if (category.labels) 
            text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                offset = offset.m[[2]], cex = text.size[[2]])
        if (!is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval1 <- zoom(zoomval)
            plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
                xlim = zoomval1[1:2], ylim = zoomval1[3:4], xlab = "", 
                ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                yaxs = "i", main = Title)
            points(Z.0, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            points(Z, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size)
            points(0, 0, pch = 3)
            if (sample.labels) 
                text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                  offset = offset.m[[1]], cex = text.size[[1]])
            if (category.labels) 
                text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                  offset = offset.m[[2]], cex = text.size[[2]])
        }
        dev.new()
        par(pty = "s", mar = parplotmar)
        cat.centroids <- LmatMinOne %*% t(Gmat) %*% Z.0[, 1:2]
        plotdata <- rbind(Z.0[, 1:2], cat.centroids[, 1:2])
        lims <- c(min(plotdata), max(plotdata))
        plot(rbind(Z.0[, 1:2], cat.centroids[, 1:2]), type = "n", 
            asp = 1, xlim = lims * exp.factor, ylim = lims * 
                exp.factor, xlab = "", ylab = "", xaxt = "n", 
            yaxt = "n", xaxs = "i", yaxs = "i", main = Title)
        points(Z.0[, 1:2], pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        points(cat.centroids, pch = pch.col.points, col = column.points.col, 
            cex = column.points.size)
        points(0, 0, pch = 3)
        if (sample.labels) 
            text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                offset = offset.m[[1]], cex = text.size[[2]])
        if (category.labels) 
            text(cat.centroids, labels = Gmat.colnames, pos = pos.m[[2]], 
                offset = offset.m[[2]], cex = text.size[[2]])
        if (!is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval1 <- zoom(zoomval)
            plot(rbind(Z.0[, 1:2], cat.centroids[, 1:2]), type = "n", 
                asp = 1, xlim = zoomval1[1:2], ylim = zoomval1[3:4], 
                xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
                xaxs = "i", yaxs = "i", main = Title)
            points(Z.0[, 1:2], pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            points(cat.centroids, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size)
            points(0, 0, pch = 3)
            if (sample.labels) 
                text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                  offset = offset.m[[1]], cex = text.size[[2]])
            if (category.labels) 
                text(cat.centroids, labels = Gmat.colnames, pos = pos.m[[2]], 
                  offset = offset.m[[2]], cex = text.size[[2]])
        }
    }
    if (mca.variant == "Burt") {
        par(pty = "s", mar = parplotmar)
        CLPs <- LmatMinHalf %*% Vmat %*% SigmaSquareMat
        CLPlist <- vector("list", p)
        varcats.names <- lapply(Glist, function(x) dimnames(x)[[2]])
        varcats.n <- lapply(varcats.names, length)
        right <- cumsum(unlist(varcats.n))
        left <- c(1, (right + 1)[-p])
        for (i in 1:p) {
            CLPlist[[i]] <- CLPs[left[i]:right[i], ]
            rownames(CLPlist[[i]]) <- dimnames(Glist[[i]])[[2]]
        }
        Z <- CLPs[, e.vects]
        Z.0 <- Gmat %*% Z/p
        plotdata <- rbind(Z[, 1:2], Z.0[, 1:2])
        lims <- c(min(plotdata), max(plotdata))
        plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
            xlim = lims * exp.factor, ylim = lims * exp.factor, 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
            yaxs = "i", main = Title)
        points(Z.0, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        points(Z, pch = pch.col.points, col = column.points.col, 
            cex = column.points.size)
        points(0, 0, pch = 3)
        if (sample.labels) 
            text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                offset = offset.m[[1]], cex = text.size[[1]])
        if (category.labels) 
            text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                offset = offset.m[[2]], cex = text.size[[2]])
        if (!is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval <- zoom(zoomval)
            plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", asp = 1, 
                xlim = zoomval[1:2], ylim = zoomval[3:4], xlab = "", 
                ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", 
                yaxs = "i", main = Title)
            points(Z.0, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            points(Z, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size)
            points(0, 0, pch = 3)
            if (sample.labels) 
                text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                  offset = offset.m[[1]], cex = text.size[[1]])
            if (category.labels) 
                text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                  offset = offset.m[[2]], cex = text.size[[2]])
        }
        if (prediction.regions) {
            dev.new()
            for (i in 1:length(CLPlist)) {
                par(pty = "s", mar = parplotmar)
                plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", 
                  asp = 1, xlim = lims * exp.factor, ylim = lims * 
                    exp.factor, xlab = "", ylab = "", xaxt = "n", 
                  yaxt = "n", xaxs = "i", yaxs = "i", main = main[i])
                pred.regions(CLPlist[[i]], x.grid = x.grid, y.grid = y.grid, 
                  plot.symbol = pred.region.symbol, size = pred.region.symbol.size, 
                  colours = region.colours)
                points(Z.0, pch = pch.row.points, col = row.points.col, 
                  cex = row.points.size)
                points(0, 0, pch = 3)
                if (sample.labels) 
                  text(Z.0, labels = dimnames(X)[[1]], pos = 1, 
                    cex = text.size[[1]])
                legend(x = "bottomleft", legend = unlist(dimnames(CLPlist[[i]])[1]), 
                  pch = rep(15, 4), col = region.colours, cex = text.size[[1]], 
                  ...)
                dev.new()
            }
        }
    }
    if (mca.variant == "EMC") {
        par(pty = "s", mar = parplotmar)
        svd.EMC <- svd((diag(n) - matrix(1, nrow = n, ncol = n)/n) %*% 
            Gmat, nu = max(n, ncat), nv = max(n, ncat))
        Sigma.EMC <- diag(svd.EMC$d)
        Vmat <- svd.EMC$v
        Umat <- svd.EMC$u
        if (ncol(Umat) > ncol(Vmat)) 
            Umat <- Umat[, 1:ncol(Vmat)]
        CLPs <- (diag(nrow(Lmat)) - matrix(1, nrow = nrow(Lmat), 
            ncol = nrow(Lmat)) %*% Lmat/n) %*% Vmat
        Z.0 <- Umat %*% Sigma.EMC
        Z <- CLPs[, e.vects]
        Z.0 <- Z.0[, e.vects]
        plotdata <- rbind(Z, Z.0)
        lims <- c(min(plotdata), max(plotdata))
        plot(rbind(Z, Z.0), type = "n", asp = 1, xlim = lims * 
            exp.factor, ylim = lims * exp.factor, xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
            main = Title)
        points(Z, pch = pch.col.points, col = column.points.col, 
            cex = column.points.size)
        points(Z.0, pch = pch.row.points, col = row.points.col, 
            cex = row.points.size)
        points(0, 0, pch = 3)
        if (sample.labels) 
            text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                offset = offset.m[[1]], cex = text.size[[1]])
        if (category.labels) 
            text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                offset = offset.m[[2]], cex = text.size[[2]])
        if (!is.null(zoomval)) {
            if (!is.numeric(zoomval)) 
                stop("zoomval must be numeric")
            zoomval <- zoom(zoomval)
            plot(rbind(Z, Z.0), type = "n", asp = 1, xlim = zoomval[1:2], 
                ylim = zoomval[3:4], xlab = "", ylab = "", xaxt = "n", 
                yaxt = "n", xaxs = "i", yaxs = "i", main = Title)
            points(Z, pch = pch.col.points, col = column.points.col, 
                cex = column.points.size)
            points(Z.0, pch = pch.row.points, col = row.points.col, 
                cex = row.points.size)
            points(0, 0, pch = 3)
            if (sample.labels) 
                text(Z.0, labels = dimnames(X)[[1]], pos = pos.m[[1]], 
                  offset = offset.m[[1]], cex = text.size[[1]])
            if (category.labels) 
                text(Z, labels = Gmat.colnames, pos = pos.m[[2]], 
                  offset = offset.m[[2]], cex = text.size[[2]])
        }
        if (prediction.regions) {
            dev.new()
            CLPlist <- vector("list", p)
            varcats.names <- lapply(Glist, function(x) dimnames(x)[[2]])
            varcats.n <- lapply(varcats.names, length)
            right <- cumsum(unlist(varcats.n))
            left <- c(1, (right + 1)[-p])
            for (i in 1:p) {
                CLPlist[[i]] <- CLPs[left[i]:right[i], ]
                rownames(CLPlist[[i]]) <- dimnames(Glist[[i]])[[2]]
            }
            for (i in 1:length(CLPlist)) {
                par(pty = "s", mar = parplotmar)
                plot(rbind(Z[, 1:2], Z.0[, 1:2]), type = "n", 
                  asp = 1, xlim = lims * exp.factor, ylim = lims * 
                    exp.factor, xlab = "", ylab = "", xaxt = "n", 
                  yaxt = "n", xaxs = "i", yaxs = "i", main = main[i])
                pred.regions(CLPlist[[i]], x.grid = x.grid, y.grid = y.grid, 
                  plot.symbol = pred.region.symbol, size = pred.region.symbol.size, 
                  colours = region.colours)
                points(Z.0, pch = pch.row.points, col = row.points.col)
                points(0, 0, pch = 3)
                if (sample.labels) 
                  text(Z.0, labels = dimnames(X)[[1]], pos = 1, 
                    cex = text.size[[1]])
                legend(x = "bottomleft", legend = unlist(dimnames(CLPlist[[i]])[1]), 
                  pch = pch.col.points, col = region.colours, 
                  cex = text.size[[1]], ...)
                dev.new()
            }
        }
    }
    if (mca.variant == "indicator" | mca.variant == "Burt") {
        e.values <- diag(SigmaSquareMat)
    }
    if (mca.variant == "EMC") {
        e.values <- diag(Sigma.EMC^2)
    }
    list(Gmat = Gmat, BurtMat = BurtMat, EMC.prop.match = EMC.prop.match, 
        Z = Z, Z.0 = Z.0, CLPlist = CLPlist, e.values = e.values)
}
