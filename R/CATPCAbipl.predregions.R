CATPCAbipl.predregions <-
function (Xcat = MCA.Table.1.data[, -1], factor.type = rep("nom", 
    ncol(Xcat)), ord.col = rep(2, 10), nom.col = 3:12, r = 2, 
    epsilon = 1e-06, orthog.transx = rep(0, ncol(Xcat)), orthog.transy = rep(0, 
        ncol(Xcat)), parplotmar = rep(3, 4), exp.factor = 1.2, 
    reverse = FALSE, select.origin = FALSE, w.factor = 2, predict.sample = NULL, 
    ort.lty = 1, offset = rep(0.5, 4), pos = c("Orthog", "Hor", 
        "Paral"), pos.text = 2, prediction.regions = "ALL") 
{
    if (!all(factor.type %in% c("nom", "ord"))) 
        stop("factor.type is a vector of 'nom' and 'ord' \n")
    n <- nrow(Xcat)
    p <- ncol(Xcat)
    pos <- pos[1]
    Gk <- vector("list", p)
    Lk.mat <- vector("list", p)
    Lk <- rep(NA, p)
    zlist <- vector("list", p)
    Hmat <- matrix(nrow = n, ncol = p)
    CatLabels <- NULL
    for (k in 1:p) {
        Gk[[k]] <- indmat(Xcat[, k])
        CatLabels <- append(CatLabels, dimnames(Gk[[k]])[[2]])
        Lk[k] <- ncol(Gk[[k]])
        Lk.mat[[k]] <- t(Gk[[k]]) %*% Gk[[k]]
    }
    G <- matrix(unlist(Gk), nrow = n)
    L.mat <- diag(diag(t(G) %*% G))
    L.mat.min.half <- diag(diag(L.mat)^-0.5)
    N.mat <- diag(n) - matrix(1, nrow = n, ncol = n)/n
    svd.mca <- svd(L.mat.min.half %*% t(G) %*% N.mat %*% G %*% 
        L.mat.min.half/p)
    z.ini <- svd.mca$v[, 1]
    names(z.ini) <- CatLabels
    count <- 0
    for (k in 1:p) {
        zlist[[k]] <- as.vector(z.ini[(count + 1):(count + ncol(Gk[[k]]))])
        names(zlist[[k]]) <- CatLabels[(count + 1):(count + ncol(Gk[[k]]))]
        count <- count + ncol(Gk[[k]])
        plot(1:length(zlist[[k]]), zlist[[k]], type = "l")
        if (identical(factor.type[k], "ord")) {
            zlist[[k]] <- ties(d = 1:length(zlist[[k]]), delta = zlist[[k]])
            lines(1:length(zlist[[k]]), zlist[[k]], col = "red")
        }
        if (identical(factor.type[k], "nom")) {
            lines(1:length(zlist[[k]]), zlist[[k]], col = "green")
        }
        dev.new()
    }
    for (k in 1:p) {
        labs <- names(zlist[[k]])
        if (identical(factor.type[k], "nom")) 
            zlist[[k]] <- as.vector(zlist[[k]]/sqrt(t(zlist[[k]]) %*% 
                t(Gk[[k]]) %*% N.mat %*% Gk[[k]] %*% zlist[[k]]))
        if (identical(factor.type[k], "ord")) {
            u.temp <- sum(diag(Lk.mat[[k]]))
            v.temp <- sum(zlist[[k]] * diag(Lk.mat[[k]]))
            w.temp <- t(zlist[[k]]) %*% Lk.mat[[k]] %*% zlist[[k]]
            b.temp <- sqrt(1/(w.temp - v.temp^2/u.temp))
            zlist[[k]] <- b.temp * zlist[[k]] - b.temp * v.temp/u.temp
        }
        names(zlist[[k]]) <- labs
        Hmat[, k] <- Gk[[k]] %*% zlist[[k]]
    }
    Hmat <- as.matrix(Hmat)
    dimnames(Hmat) <- dimnames(Xcat)
    again <- TRUE
    RSS.old <- NULL
    Hmat.cent <- N.mat %*% Hmat
    iter <- 0
    while (again) {
        iter <- iter + 1
        svd.out <- svd(Hmat.cent)
        Yr <- scale(svd.out$u[, 1:r] %*% diag(svd.out$d[1:r]) %*% 
            t(svd.out$v[, 1:r]), scale = FALSE, center = TRUE)
        RSS <- sum((Hmat.cent - Yr)^2)
        if (!is.null(RSS.old)) 
            again <- ((RSS.old - RSS) > epsilon)
        RSS.old <- RSS
        if (again) {
            for (k in 1:p) {
                labs <- names(zlist[[k]])
                if (identical(factor.type[k], "nom")) {
                  scale.factor <- as.numeric(sqrt(t(Yr[, k]) %*% 
                    Gk[[k]] %*% solve(Lk.mat[[k]]) %*% t(Gk[[k]]) %*% 
                    Yr[, k]))
                  zlist[[k]] <- as.vector((solve(Lk.mat[[k]]) %*% 
                    t(Gk[[k]]) %*% Yr[, k])/scale.factor)
                }
                if (identical(factor.type[k], "ord")) {
                  zlist[[k]] <- as.vector((solve(Lk.mat[[k]]) %*% 
                    t(Gk[[k]]) %*% Yr[, k]))
                  zlist[[k]] <- ties(d = 1:length(zlist[[k]]), 
                    delta = zlist[[k]])
                  u.temp <- sum(diag(Lk.mat[[k]]))
                  v.temp <- sum(zlist[[k]] * diag(Lk.mat[[k]]))
                  w.temp <- t(zlist[[k]]) %*% Lk.mat[[k]] %*% 
                    zlist[[k]]
                  b.temp <- sqrt(1/(w.temp - v.temp^2/u.temp))
                  zlist[[k]] <- b.temp * zlist[[k]] - b.temp * 
                    v.temp/u.temp
                }
                names(zlist[[k]]) <- labs
                Hmat[, k] <- Gk[[k]] %*% zlist[[k]]
            }
        }
        Hmat.cent <- N.mat %*% Hmat
    }
    Lz <- rep(NA, p)
    for (i in 1:p) Lz[i] <- sum(Lk.mat[[i]] %*% zlist[[i]])
    zLz <- rep(NA, p)
    for (i in 1:p) zLz[i] <- sum(t(zlist[[i]]) %*% Lk.mat[[i]] %*% 
        zlist[[i]])
    Hmat.cent.svd <- svd(Hmat.cent)
    par(pty = "s", mar = parplotmar)
    vars.points <- Hmat.cent.svd$v
    samples.points <- Hmat.cent %*% Hmat.cent.svd$v[, 1:2]
    grad <- rep(NA, p)
    for (i in 1:p) {
        regr.line <- lm(c(0, vars.points[i, 2]) ~ c(0, vars.points[i, 
            1]))
        grad[i] <- coefficients(regr.line)[2]
    }
    means <- apply(Hmat.cent, 2, mean)
    sd <- sqrt(apply(Hmat.cent, 2, var))
    axes.rows <- 1/(diag(vars.points %*% t(vars.points))) * vars.points
    z.axes <- lapply(1:p, function(j, X.1, axes.rows, orthog.transx, 
        orthog.transy) {
        phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
            axes.rows %*% c(orthog.transx[j], orthog.transy[j], 
            rep(0, p - 2))
        axis.vals <- (X.1[, j])
        axis.points <- matrix(0, nrow = length(axis.vals), ncol = p + 
            2)
        axis.points[, 1] <- orthog.transx[j] + (axis.vals - phi.vec[j]) * 
            axes.rows[j, 1]
        axis.points[, 2] <- orthog.transy[j] + (axis.vals - phi.vec[j]) * 
            axes.rows[j, 2]
        for (k in 3:p) axis.points[, k] <- axis.vals * axes.rows[j, 
            k]
        axis.points[, p + 1] <- axis.vals
        axis.points[, p + 2] <- NA
        axis.points <- axis.points[order(axis.points[, 1]), ]
        return(unique(zapsmall(axis.points)))
    }, X.1 = Hmat.cent, axes.rows = axes.rows, orthog.transx = orthog.transx, 
        orthog.transy = orthog.transy)
    all.axes.points <- rbind(z.axes[[1]][, 1:2], z.axes[[2]][, 
        1:2])
    for (i in 3:p) all.axes.points <- rbind(all.axes.points, 
        z.axes[[i]][, 1:2])
    all.points <- rbind(all.axes.points, samples.points)
    plot(all.points[, 1] * exp.factor, all.points[, 2] * exp.factor, 
        xlim = range(all.points[, 1] * exp.factor), ylim = range(all.points[, 
            2] * exp.factor), xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", type = "n", xaxs = "i", yaxs = "i", asp = 1)
    for (i in 1:p) {
        z.axes[[i]] <- as.data.frame(z.axes[[i]])
        if (factor.type[i] == "nom") 
            z.axes[[i]][, p + 2] <- names(zlist[[i]])[match(round(z.axes[[i]][, 
                p + 1], digits = 7), round(zlist[[i]], digits = 7))]
        if (factor.type[i] == "ord") {
            labs.temp <- tied.names(names(zlist[[i]]), zlist[[i]])
            zlist.temp <- unique(zlist[[i]])
            names(zlist.temp) <- labs.temp
            z.axes[[i]][, p + 2] <- names(zlist.temp)[match(round(z.axes[[i]][, 
                p + 1], digits = 7), round(zlist.temp, digits = 7))]
            labs.temp <- NULL
            Zlist.temp <- NULL
        }
    }
    if (!is.null(prediction.regions)) {
        for (i in 1:p) {
            plot(all.points[, 1] * exp.factor, all.points[, 2] * 
                exp.factor, xlim = range(all.points[, 1] * exp.factor), 
                ylim = range(all.points[, 2] * exp.factor), xaxt = "n", 
                yaxt = "n", xlab = "", ylab = "", type = "n", 
                xaxs = "i", yaxs = "i", asp = 1)
            pred.regions(z.axes[[i]][, 1:p], 0.005, 0.005, 0.5)
            dev.new()
        }
    }
    for (i in 1:p) {
        Draw.line2(x.vals = z.axes[[i]][, 1], y.vals = z.axes[[i]][, 
            2], marker.vals = z.axes[[i]][, 1], line.name = dimnames(Xcat)[[2]][i], 
            axis.col = "white", pos = pos, offset = offset)
        n.cat <- nrow(z.axes[[i]])
        if (factor.type[i] == "nom") 
            draw.halfway.lines(x = z.axes[[i]][, 1], y = z.axes[[i]][, 
                2], colour = nom.col[1:n.cat], width = rep(w.factor, 
                n.cat))
        if (factor.type[i] == "ord") {
            temp <- z.axes[[i]][, 1:2]
            draw.halfway.lines(x = temp[, 1], y = temp[, 2], 
                colour = ord.col[1:n], width = 1:(n.cat + 1), 
                w.factor = w.factor, reverse = reverse)
        }
        points(z.axes[[i]][, 1:2], pch = 16)
        for (j in 1:n.cat) {
            text(x = z.axes[[i]][j, 1], y = z.axes[[i]][j, 2], 
                label = z.axes[[i]][j, p + 2], cex = 0.65, col = "red", 
                pos = pos.text, offset = 0.4)
        }
    }
    points(samples.points, pch = 15, col = "black")
    text(samples.points, labels = dimnames(Xcat)[[1]], pos = pos.text, 
        cex = 0.65, offset = 0.4)
    if (!is.null(predict.sample)) 
        for (i in 1:p) DrawOrthogline(x1 = z.axes[[i]][1, 1], 
            x2 = z.axes[[i]][2, 1], y1 = z.axes[[i]][1, 2], y2 = z.axes[[i]][2, 
                2], px = samples.points[predict.sample, 1], py = samples.points[predict.sample, 
                2], ort.lty = ort.lty)
    if (select.origin) {
        cat("Move cursor where you want origin and press left mouse button \n")
        flush.console()
        origin.pos <- locator(1)
        dev.off()
        z.axes <- lapply(1:p, function(j, X.1, axes.rows, orthog.transx, 
            orthog.transy) {
            phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                axes.rows %*% c(orthog.transx[j], orthog.transy[j], 
                rep(0, p - 2))
            axis.vals <- (X.1[, j])
            axis.points <- matrix(0, nrow = length(axis.vals), 
                ncol = p + 2)
            axis.points[, 1] <- orthog.transx[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 1]
            axis.points[, 2] <- orthog.transy[j] + (axis.vals - 
                phi.vec[j]) * axes.rows[j, 2]
            for (k in 3:p) axis.points[, k] <- axis.vals * axes.rows[j, 
                k]
            axis.points[, p + 1] <- axis.vals
            axis.points[, p + 2] <- NA
            axis.points <- axis.points[order(axis.points[, 1]), 
                ]
            return(unique(zapsmall(axis.points)))
        }, X.1 = Hmat.cent, axes.rows = axes.rows, orthog.transx = rep(origin.pos$x, 
            p), orthog.transy = rep(origin.pos$y, p))
        all.axes.points <- rbind(z.axes[[1]][, 1:2], z.axes[[2]][, 
            1:2])
        for (i in 3:p) all.axes.points <- rbind(all.axes.points, 
            z.axes[[i]][, 1:2])
        all.points <- rbind(all.axes.points, samples.points)
        plot(all.points[, 1] * exp.factor, all.points[, 2] * 
            exp.factor, xlim = range(all.points[, 1] * exp.factor), 
            ylim = range(all.points[, 2] * exp.factor), xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", type = "n", xaxs = "i", 
            yaxs = "i", asp = 1)
        for (i in 1:p) {
            z.axes[[i]] <- as.data.frame(z.axes[[i]])
            if (factor.type[i] == "nom") 
                z.axes[[i]][, p + 2] <- names(zlist[[i]])[match(round(z.axes[[i]][, 
                  p + 1], digits = 7), round(zlist[[i]], digits = 7))]
            if (factor.type[i] == "ord") {
                labs.temp <- tied.names(names(zlist[[i]]), zlist[[i]])
                zlist.temp <- unique(zlist[[i]])
                names(zlist.temp) <- labs.temp
                z.axes[[i]][, p + 2] <- names(zlist.temp)[match(round(z.axes[[i]][, 
                  p + 1], digits = 7), round(zlist.temp, digits = 7))]
                labs.temp <- NULL
                Zlist.temp <- NULL
            }
        }
        for (i in 1:p) {
            Draw.line2(x.vals = z.axes[[i]][, 1], y.vals = z.axes[[i]][, 
                2], marker.vals = z.axes[[i]][, p + 1], line.name = dimnames(Xcat)[[2]][i], 
                axis.col = "white", pos = pos, offset = offset)
            n.cat <- nrow(z.axes[[i]])
            if (factor.type[i] == "nom") 
                draw.halfway.lines(x = z.axes[[i]][, 1], y = z.axes[[i]][, 
                  2], colour = nom.col[1:n.cat], width = rep(w.factor, 
                  n.cat))
            if (factor.type[i] == "ord") {
                temp <- z.axes[[i]][, 1:2]
                draw.halfway.lines(x = temp[, 1], y = temp[, 
                  2], colour = ord.col[1:n], width = 1:(n.cat + 
                  1), w.factor = w.factor, reverse = reverse)
            }
            points(z.axes[[i]][, 1:2], pch = 16)
            for (j in 1:n.cat) {
                text(x = z.axes[[i]][j, 1], y = z.axes[[i]][j, 
                  2], label = z.axes[[i]][j, p + 2], cex = 0.65, 
                  col = "red", pos = pos.text, offset = 0.4)
            }
        }
        points(samples.points, pch = 15, col = "black")
        text(samples.points, labels = dimnames(Xcat)[[1]], pos = pos.text, 
            cex = 0.65, offset = 0.4)
        if (!is.null(predict.sample)) 
            for (i in 1:p) DrawOrthogline(x1 = z.axes[[i]][1, 
                1], x2 = z.axes[[i]][2, 1], y1 = z.axes[[i]][1, 
                2], y2 = z.axes[[i]][2, 2], px = samples.points[predict.sample, 
                1], py = samples.points[predict.sample, 2], ort.lty = ort.lty)
    }
    dev.new()
    PCAbipl.cat(Hmat.cent, exp.factor = exp.factor)
    for (k in 1:p) {
        dev.new()
        plot(1:length(zlist[[k]]), zlist[[k]], type = "b", xaxt = "n", 
            ylab = "Final z quantification", xlab = dimnames(Xcat)[[2]][k])
        axis(side = 1, at = 1:length(names(zlist[[k]])), labels = names(zlist[[k]]))
    }
    list(iter = iter, Hmat = Hmat, Hmat.cent = Hmat.cent, Yr = Yr, 
        zlist = zlist, z.axes = z.axes)
}
