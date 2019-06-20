CVAbipl.3d <-
function (X = Ocotea.data[, 3:8], G = indmat(Ocotea.data[, 2]), 
    X.new.samples = NULL, e.vects = 1:3, weightedCVA = c("weighted", 
        "unweightI", "unweightCent"), ax.type = c("predictive", 
        "interpolative"), alpha.3d = 0.7, ax.col = rep("black", 
        ncol(X)), ax.col.3d = "black", aspect.3d = "iso", cex.3d = 0.6, 
    colour.scheme = NULL, colours = c(4:12, 3:1), colours.3d.means = c(4:12, 
        3:1), col.plane.3d = "lightgrey", col.text.3d = "black", 
    factor.x = 2, factor.y = 2, factor.3d.axes = 0.01, font.3d = 2, 
    ID.labs = FALSE, ID.3d = 1:nrow(X), n.int = rep(5, ncol(X)), 
    predictions.sample = NULL, size.ax.3d = 0.5, size.points.3d = 0.2, 
    size.means.3d = 5, specify.classes = dimnames(G)[[2]], Titles.3d = c("", 
        "", "x", "y", "z")) 
{
    dim.biplot <- 3
    ax.type <- ax.type[1]
    weightedCVA <- weightedCVA[1]
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    if (is.null(G)) {
        cat("Use function PCAbipl for a PCA biplot or provide G for a CVA biplot\n")
        return(invisible())
    }
    if (ncol(G) < 4) 
        warning("Three class CVA biplot reduces to a plane;  third axis not uniquely determined \n\n             Two class CVA biplot reduces to a line;  second and third axes not uniquely determined \n")
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    if (is.null(dimnames(G))) 
        dimnames(G) <- list(NULL, paste("class", 1:J, sep = ""))
    if (length(dimnames(G)[[2]]) == 0) 
        dimnames(G)[[2]] <- paste("class", 1:J, sep = "")
    n.groups <- apply(G, 2, sum)
    if (is.numeric(specify.classes)) 
        specify.classes <- dimnames(G)[[2]][specify.classes]
    if (!is.null(specify.classes)) {
        if (!is.numeric(specify.classes)) {
            test <- match(specify.classes, unique(as.character(dimnames(G)[[2]])), 
                nomatch = 0)
            if (sum(test != 0) < length(specify.classes)) 
                stop(" One or more of specified classes non-existent.  Check spelling")
        }
    }
    if (!is.null(specify.classes)) {
        select.numeric.classes <- rep(0, length(specify.classes))
        for (k in 1:length(specify.classes)) select.numeric.classes[k] <- (1:length(dimnames(G)[[2]]))[dimnames(G)[[2]] == 
            specify.classes[k]]
    }
    else {
        if (is.null(specify.classes)) 
            select.numeric.classes <- NULL
    }
    dimnames(G)[[2]] <- paste(dimnames(G)[[2]], "; n = ", n.groups, 
        sep = "")
    X <- as.matrix(X)
    unscaled.X <- X
    X.cent <- scale(X, center = TRUE, scale = FALSE)
    X.means <- apply(X, 2, mean)
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    K <- min(p, J - 1)
    if (nrow(X) != nrow(G)) {
        stop("Rows of X and G unequal!!")
    }
    class.vec <- apply(t(apply(G, 1, function(x) x == max(x))), 
        1, function(s, G) dimnames(G)[[2]][s], G = G)
    if (J > length(colours)) 
        stop(paste("Increase size of colours argument"))
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste("s", 1:n, sep = ""), paste("V", 
            1:p, sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste("s", 1:n, sep = "")
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
    Gmat <- G
    Nmat <- t(Gmat) %*% Gmat
    XcentBar.groups <- solve(Nmat) %*% t(Gmat) %*% X.cent
    XBar.groups <- solve(Nmat) %*% t(Gmat) %*% X
    SSP.T <- t(X.cent) %*% X.cent
    SSP.B <- t(XcentBar.groups) %*% Nmat %*% XcentBar.groups
    SSP.W <- SSP.T - SSP.B
    Wmat <- SSP.W
    svd.Wmat <- svd(Wmat)
    lambdamatI <- diag(svd.Wmat$d)
    Lmat <- svd.Wmat$u %*% solve(sqrt(lambdamatI))
    if (weightedCVA == "weighted") {
        Cmat <- Nmat
        Csqrt <- sqrt(Nmat)
    }
    if (weightedCVA == "unweightedI") {
        Cmat <- diag(J)
        Csqrt <- Cmat
    }
    if (weightedCVA == "unweightedCent") {
        Cmat <- diag(J) - matrix(1/J, nrow = J, ncol = J)
        Csqrt <- Cmat
    }
    if (is.na(match(weightedCVA, c("weighted", "unweightedI", 
        "unweightedCent")))) 
        stop(" Argument 'weightedCVA' must be one of 'weighted','unweightedI','unweightedCent' ")
    svd.step2 <- svd(t(Lmat) %*% t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups %*% Lmat)
    Vmat <- svd.step2$v
    lambdamat <- diag(svd.step2$d)
    svd.2sided <- Eigen.twosided(t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups, Wmat)
    Mmat <- svd.2sided$W
    lambdamat.2sided <- svd.2sided$Lambda.mat
    vec.temp <- rep(0, p)
    vec.temp[1:dim.biplot] <- 1
    Jmat <- diag(vec.temp)
    XcentBarLVJ <- XcentBar.groups %*% Mmat %*% Jmat
    XcentLVJ <- X.cent %*% Mmat %*% Jmat
    XcentBarHat <- XcentBar.groups %*% Mmat %*% Jmat %*% solve(Mmat)
    if (length(e.vects) == 3) 
        CVA.mean <- XcentBar.groups %*% Mmat[, e.vects[1:3]]
    else stop("e.vects must be specified to have 3 elements. /n")
    Z.means.mat <- data.frame(CVA.mean, colr = colours.3d.means[1:J], 
        stringsAsFactors = FALSE)
    classnames <- dimnames(Z.means.mat)[[1]]
    Z <- X.cent %*% Mmat[, e.vects[1:3]]
    Z <- data.frame(Z, colr = as.character(colours[1]), stringsAsFactors = FALSE)
    if (J > 0) 
        for (j in 1:J) {
            Z[G[, j] == 1, 4] <- as.character(colours[j])
        }
    Brr <- solve(Mmat)[e.vects[1:3], ]
    Br <- (Mmat)[, e.vects[1:3]]
    if (ax.type == "predictive") 
        axes.rows <- solve(diag(diag(t(Brr) %*% Brr))) %*% t(Brr)
    else {
        if (ax.type == "interpolative") 
            axes.rows <- Br
        else stop("ax.type is either 'interpolative' or 'predictive' (the default). \n")
    }
    axes.rows <- solve(diag(diag(t(Brr) %*% Brr))) %*% t(Brr)
    z.axes <- lapply(1:p, function(j, unscaled.X, means, axes.rows, 
        n.int) {
        number.points <- 50
        std.markers <- pretty(unscaled.X[, j], n = n.int[j])
        std.range <- c(min(std.markers), max(std.markers))
        std.markers.min <- std.markers - (std.range[2] - std.range[1]) * 
            factor.3d.axes
        std.markers.max <- std.markers + (std.range[2] - std.range[1]) * 
            factor.3d.axes
        std.markers <- c(std.markers, std.markers.min, std.markers.max)
        interval <- std.markers - means[j]
        axis.vals <- seq(from = min(interval), to = max(interval), 
            length = number.points)
        axis.vals <- sort(unique(c(axis.vals, interval)))
        number.points <- length(axis.vals)
        axis.points <- matrix(0, nrow = number.points, ncol = 5)
        axis.points[, 1] <- axis.vals * axes.rows[j, 1]
        axis.points[, 2] <- axis.vals * axes.rows[j, 2]
        axis.points[, 3] <- axis.vals * axes.rows[j, 3]
        axis.points[, 4] <- axis.vals + means[j]
        axis.points[, 5] <- 0
        for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
            4] - std.markers) == 0)) 
            axis.points[i, 5] <- 1
        return(axis.points)
    }, unscaled.X = unscaled.X, means = X.means, axes.rows = axes.rows, 
        n.int = n.int)
    AxisPredictivity <- diag(t(XcentBarHat) %*% Cmat %*% XcentBarHat)/diag(t(XcentBar.groups) %*% 
        Cmat %*% XcentBar.groups)
    ClassPredictivity <- diag(Csqrt %*% XcentBarHat %*% solve(Wmat) %*% 
        t(XcentBarHat) %*% Csqrt)/diag(Csqrt %*% XcentBar.groups %*% 
        solve(Wmat) %*% t(XcentBar.groups) %*% Csqrt)
    X.within <- (diag(n) - Gmat %*% (solve(Nmat)) %*% t(Gmat)) %*% 
        X.cent
    X.within.hat <- X.within %*% Mmat %*% Jmat %*% solve(Mmat)
    Within.group.axis.predictivity <- diag(t(X.within.hat) %*% 
        X.within.hat)/(diag(t(X.within) %*% X.within))
    Within.group.sample.predictivity <- diag((X.within.hat) %*% 
        solve(Wmat) %*% t(X.within.hat))/(diag((X.within) %*% 
        solve(Wmat) %*% t(X.within)))
    if (is.null(G)) 
        stop("G cannot be set to NULL. Groups must be specified for a CVA biplot \n")
    if (!is.null(colour.scheme)) {
        my.colours <- colorRampPalette(colour.scheme)
        colours <- my.colours(colours)
    }
    lambda.vec <- zapsmall(diag(lambdamat))
    CVA.quality.canvar <- sum(lambda.vec[e.vects[1:dim.biplot]])/sum(lambda.vec)
    CVA.quality.origvar <- sum(diag(t(XcentBarHat) %*% Cmat %*% 
        XcentBarHat))/sum(diag(t(XcentBar.groups) %*% Cmat %*% 
        XcentBar.groups))
    adequacy.p <- diag(Br %*% t(Br))/diag(Mmat %*% t(Mmat))
    names(adequacy.p) <- dimnames(X)[[2]]
    {
        require(rgl)
        z.axes.names <- dimnames(X)[[2]]
        graph.coord <- z.axes[[1]][c(1, nrow(z.axes[[1]])), 1:3]
        specify.classes <- select.numeric.classes
        plot(1:10, 1:10, type = "n", xlab = "", ylab = "", yaxt = "n", 
            xaxt = "n", bty = "n")
        legend(x = 5, y = 5, legend = dimnames(Z.means.mat)[[1]], 
            pch = 15, col = Z.means.mat[, 4], cex = 1.5)
        .samples.plot.3d <- function(Z, Z.means.mat, class.vec, 
            specify.classes, size.points.3d) {
            if (is.null(Z.means.mat)) 
                J <- 0
            else J <- nrow(Z.means.mat)
            classes <- unique(dimnames(Z.means.mat)[[1]])[specify.classes]
            for (j in classes) {
                Z.class <- data.frame(Z[class.vec == j, , drop = FALSE], 
                  stringsAsFactors = FALSE)
                spheres3d(Z.class[, 1:3], radius = size.points.3d, 
                  color = Z.class[, 4])
            }
        }
        for (i in 2:p) graph.coord <- rbind(graph.coord, z.axes[[i]][c(1, 
            nrow(z.axes[[i]])), 1:3])
        open3d()
        text3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
        points3d(Z[, 1:3], size = 0)
        aspect3d(aspect.3d)
        if (!is.null(specify.classes)) 
            .samples.plot.3d(Z = Z, Z.means.mat = Z.means.mat, 
                class.vec = class.vec, specify.classes = specify.classes, 
                size.points.3d = size.points.3d)
        if (ID.labs) 
            text3d(Z[, 1:3], text = ID.3d, col = "black", font = font.3d, 
                cex = cex.3d)
        for (i in 1:p) {
            axis <- z.axes[[i]]
            text3d(axis[axis[, 4] == max(axis[, 4]), 1:3], text = z.axes.names[i], 
                col = col.text.3d, font = font.3d, cex = cex.3d)
            axis <- axis[axis[, 4] != max(axis[, 4]), 1:3]
            lines3d(axis[c(1, nrow(axis)), 1:3], size = size.ax.3d, 
                color = ax.col[i])
        }
        points3d(Z.means.mat[, 1:3], size = size.means.3d, color = Z.means.mat[, 
            4])
        if (!(is.null(X.new.samples))) {
            pch.new <- paste("N", 1:nrow(X.new.samples), sep = "")
            Z.new <- scale(X.new.samples, means, scale = FALSE) %*% 
                Br
            text3d(Z.new, text = pch.new, adj = c(0.5, 0.5), 
                col = "black", font = font.3d, cex = cex.3d)
        }
        surface3d(x = c(min(graph.coord[, 1]), max(graph.coord[, 
            1])) * factor.x, y = c(min(graph.coord[, 2]), max(graph.coord[, 
            2])) * factor.y, z = matrix(0, nrow = 2, ncol = 2), 
            color = col.plane.3d, alpha = alpha.3d)
        axes3d(col = ax.col.3d)
        title3d(main = Titles.3d[1], sub = Titles.3d[2], xlab = Titles.3d[3], 
            ylab = Titles.3d[4], zlab = Titles.3d[5])
    }
    device.ID <- rgl.cur()
    answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    while (!(answer == "Y" | answer == "y" | answer == "N" | 
        answer == "n")) answer <- readline("Save 3D graph as a .png file? Y/N  \n")
    if (answer == "Y" | answer == "y") 
        repeat {
            file.name <- readline("Provide file name including full NOT in quotes and SINGLE back slashes!  \n")
            file.name <- paste(file.name, ".png", sep = "")
            snapshot3d(file = file.name)
            rgl.set(device.ID)
            answer2 <- readline("Save another 3D graph as a .png file? Y/N  \n")
            if (answer2 == "Y" | answer2 == "y") 
                next
            else break
        }
    else rgl.set(device.ID)
    centered.vals <- attr(scale(unscaled.X, scale = FALSE, center = TRUE), 
        "scaled:center")
    means.predictions.3D <- CVA.predictions.mat.3d(X = XcentBar.groups, 
        B = Mmat)$predictions.for.originalX[[3]]
    means.predictions.3D <- scale(means.predictions.3D, scale = FALSE, 
        center = -centered.vals)
    samples.predictions.3D <- CVA.predictions.mat.3d(X = unscaled.X, 
        B = Mmat)$predictions.for.originalX[[3]]
    if (is.null(predictions.sample)) 
        samples.predictions.3D <- NULL
    list(XcentBar.groups = XcentBar.groups, XBar.groups = XBar.groups, 
        Mmat = Mmat, lambda = diag(zapsmall(lambda.vec)), CVA.quality.origvar = CVA.quality.origvar, 
        CVA.quality.canvar = CVA.quality.canvar, adequacy = adequacy.p, 
        means.predictions.3D = means.predictions.3D, samples.predictions.3D = samples.predictions.3D, 
        Z.means.mat = Z.means.mat)
}
