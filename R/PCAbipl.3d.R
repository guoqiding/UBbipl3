PCAbipl.3d <-
function (X = Ocotea.data[, 3:8], G = NULL, X.new = NULL, scaled.mat = TRUE, 
    e.vects = 1:3, alpha.3d = 0.7, ax.col = rep("black", ncol(X)), 
    ax.col.3d = "black", aspect.3d = "iso", cex.3d = 0.6, ax.type = c("predictive", 
        "interpolative"), adjust.3d = c(0.5, 0.5), colour.scheme = NULL, 
    colours = c(4:12, 3:1), col.plane.3d = "lightgrey", col.text.3d = "black", 
    correlation.biplot = FALSE, factor.x = 2, factor.y = 2, font.3d = 2, 
    ID.labs = FALSE, ID.3d = rownames(X), n.int = rep(5, ncol(X)), 
    plot.class.means = FALSE, predictions = TRUE, size.ax.3d = 1, 
    size.means.3d = 10, size.points.3d = 10, specify.classes = dimnames(G)[[2]], 
    Title = NULL, Titles.3d = c("", "", "Dim 1", "Dim 2", "Dim 3")) 
{
    if (is.null(G)) {
        G <- matrix(indmat(rep(1, nrow(X))), ncol = 1)
        dimnames(G) <- list(1:nrow(X), "AllData")
    }
    ax.type <- ax.type[1]
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    if (!is.null(G)) 
        if (length(specify.classes) > ncol(G)) 
            return(cat("Number specified classes must not be larger than the number of different classes\n"))
    X <- as.matrix(X)
    unscaled.X <- X
    means <- apply(X, 2, mean)
    sd <- sqrt(apply(X, 2, var))
    if (scaled.mat) 
        X <- scale(X)
    else {
        X <- scale(X, scale = FALSE)
        sd <- rep(1, ncol(X))
    }
    n <- nrow(X)
    p <- ncol(X)
    J <- ncol(G)
    if (J > length(colours)) 
        stop(paste("Increase size of colours argument"))
    if (is.null(dimnames(X))) 
        dimnames(X) <- list(paste("s", 1:n), paste("V", 1:p, 
            sep = ""))
    if (length(dimnames(X)[[1]]) == 0) 
        dimnames(X)[[1]] <- paste("s", 1:n)
    if (length(dimnames(X)[[2]]) == 0) 
        dimnames(X)[[2]] <- paste("V", 1:p, sep = "")
    if (nrow(G) != n) 
        stop("number of rows of X and G differ")
    if (is.null(dimnames(G))) 
        dimnames(G) <- list(NULL, paste("class", 1:J, sep = ""))
    if (length(dimnames(G)[[2]]) == 0) 
        dimnames(G)[[2]] <- paste("class", 1:J, sep = "")
    n.groups <- apply(G, 2, sum)
    if (is.numeric(specify.classes)) 
        stop("Class names must be supplied, not numeric values")
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
    if (length(dimnames(G)[[2]]) == 1) 
        class.vec <- rep(dimnames(G)[[2]], nrow(X))
    else {
        class.vec <- apply(t(apply(G, 1, function(x) x == max(x))), 
            1, function(s, G) dimnames(G)[[2]][s], G = G)
    }
    svd.out <- svd(t(X) %*% X)
    V.mat <- svd.out$u
    Vr <- svd.out$u[, e.vects[1:3]]
    eigval <- svd.out$d
    lambda.mat <- diag(eigval)
    eigval.r <- eigval[e.vects[1:3]]
    lambda.r.mat <- diag(eigval.r)
    predictivity.mat <- diag(diag(Vr %*% lambda.r.mat %*% t(Vr))) %*% 
        solve(diag(diag(V.mat %*% lambda.mat %*% t(V.mat))))
    predictivities <- round(diag(predictivity.mat), digits = 3)
    names(predictivities) <- dimnames(X)[[2]]
    PCA.quality <- paste("Quality of display =", round(((eigval[e.vects[1]] + 
        eigval[e.vects[2]] + eigval[e.vects[3]])/sum(eigval)) * 
        100, digits = 2), "%")
    fit.adequacy <- round(diag(Vr %*% t(Vr)), digits = 3)
    names(fit.adequacy) <- dimnames(X)[[2]]
    {
        class.means.mat <- solve(t(G) %*% G) %*% t(G) %*% unscaled.X
    }
    if (correlation.biplot) {
        lambda.r <- diag(svd(t(X) %*% X)$d[1:3])
        Z <- sqrt(n - 1) * X %*% Vr %*% (sqrt(solve(lambda.r)))
        Z.means.mat <- sqrt(n - 1) * scale(class.means.mat, means, 
            sd) %*% Vr %*% (sqrt(solve(lambda.r)))
        Z.means.mat <- data.frame(Z.means.mat, colr = as.character(colours[1:J]), 
            stringsAsFactors = FALSE)
        classnames <- dimnames(Z.means.mat)[[1]]
    }
    else {
        Z <- X %*% Vr
        Z.means.mat <- scale(class.means.mat, means, sd) %*% 
            Vr
        Z.means.mat <- data.frame(Z.means.mat, colr = as.character(colours[1:J]), 
            stringsAsFactors = FALSE)
        classnames <- dimnames(Z.means.mat)[[1]]
    }
    Z <- data.frame(Z, colr = as.character(colours[1]), stringsAsFactors = FALSE)
    if (J > 0) 
        for (j in 1:J) Z[G[, j] == 1, 4] <- colours[j]
    if (ax.type == "predictive") 
        axes.rows <- 1/(diag(Vr %*% t(Vr))) * Vr
    else {
        if (ax.type == "interpolative") 
            axes.rows <- Vr
        else stop("ax.type is either 'interpolative' or 'predictive' (the default)")
    }
    if (correlation.biplot) {
        axes.rows <- (sqrt(n - 1)/(diag(Vr %*% lambda.r %*% t(Vr)))) * 
            Vr %*% sqrt(lambda.r)
    }
    z.axes <- lapply(1:p, function(j, unscaled.X, means, sd, 
        axes.rows, n.int) {
        number.points <- 100
        std.markers <- pretty(unscaled.X[, j], n = n.int[j])
        std.range <- c(min(std.markers), max(std.markers))
        std.markers.min <- std.markers - (std.range[2] - std.range[1])
        std.markers.max <- std.markers + (std.range[2] - std.range[1])
        std.markers <- c(std.markers, std.markers.min, std.markers.max)
        interval <- (std.markers - means[j])/sd[j]
        axis.vals <- seq(from = min(interval), to = max(interval), 
            length = number.points)
        axis.vals <- sort(unique(c(axis.vals, interval)))
        number.points <- length(axis.vals)
        axis.points <- matrix(0, nrow = number.points, ncol = 5)
        axis.points[, 1] <- axis.vals * axes.rows[j, 1]
        axis.points[, 2] <- axis.vals * axes.rows[j, 2]
        axis.points[, 3] <- axis.vals * axes.rows[j, 3]
        axis.points[, 4] <- axis.vals * sd[j] + means[j]
        axis.points[, 5] <- 0
        for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
            4] - std.markers) == 0)) 
            axis.points[i, 5] <- 1
        return(axis.points)
    }, unscaled.X = unscaled.X, means = means, sd = sd, axes.rows = axes.rows, 
        n.int = n.int)
    {
        require(rgl)
        z.axes.names = colnames(X)
        graph.coord <- z.axes[[1]][c(1, nrow(z.axes[[1]])), 1:3]
        specify.classes <- select.numeric.classes
        plot(1:10, 1:10, type = "n", xlab = "", ylab = "", yaxt = "n", 
            xaxt = "n", bty = "n")
        legend(x = 5, y = 5, legend = dimnames(Z.means.mat)[[1]], 
            pch = 15, col = Z.means.mat[, 4], cex = 1.5)
        .samples.plot.3d <- function(Z, Z.means.mat, adjust.3d, 
            class.vec, specify.classes, size.points.3d) {
            if (is.null(Z.means.mat)) 
                J <- 0
            else J <- nrow(Z.means.mat)
            classes <- unique(dimnames(Z.means.mat)[[1]])[specify.classes]
            for (j in classes) {
                Z.class <- data.frame(Z[class.vec == j, ], stringsAsFactors = FALSE)
                spheres3d(Z.class[, 1:3], radius = size.points.3d, 
                  color = Z.class[, 4])
                if (ID.labs) 
                  texts3d(Z.class[, 1:3], adj = adjust.3d, texts = ID.3d, 
                    col = Z.class[, 4], font = font.3d, cex = cex.3d)
            }
        }
        for (i in 2:p) graph.coord <- rbind(graph.coord, z.axes[[i]][c(1, 
            nrow(z.axes[[i]])), 1:3])
        open3d()
        texts3d(0, 0, 0, text = "", font = font.3d, cex = cex.3d)
        points3d(graph.coord, size = 0)
        aspect3d(aspect.3d)
        if (!is.null(specify.classes)) 
            .samples.plot.3d(Z = Z, Z.means.mat = Z.means.mat, 
                adjust.3d = adjust.3d, class.vec = class.vec, 
                specify.classes = specify.classes, size.points.3d = size.points.3d)
        for (i in 1:p) {
            axis <- z.axes[[i]]
            texts3d(axis[axis[, 4] == max(axis[, 4]), 1:3], text = z.axes.names[i], 
                col = col.text.3d, font = font.3d, cex = cex.3d)
            texts3d(axis[axis[, 4] == min(axis[, 4]), 1:3], text = z.axes.names[i], 
                col = col.text.3d, font = font.3d, cex = cex.3d)
            axis <- axis[axis[, 4] != max(axis[, 4]), 1:3]
            lines3d(axis[c(1, nrow(axis)), 1:3], size = size.ax.3d, 
                color = ax.col[i])
        }
        if (plot.class.means) 
            points3d(x = Z.means.mat[, 1], y = Z.means.mat[, 
                2], z = Z.means.mat[, 3], size = size.means.3d, 
                color = Z.means.mat[, 4])
        if (!(is.null(X.new))) {
            pch.new <- paste("N", 1:nrow(X.new), sep = "")
            Z.new <- scale(X.new, means, scale = FALSE) %*% Vr
            texts3d(Z.new, text = pch.new, adj = c(0.5, 0.5), 
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
    if (predictions) 
        predictions <- PCA.predictions.mat(X = unscaled.X, scaled.mat = scaled.mat, 
            e.vects = e.vects[1:3])
    else predictions = NULL
    list(Z = Z, z.axes = z.axes, V = V.mat, Eigenvectors = Vr, 
        e.vals = eigval, PCA.quality = PCA.quality, adequacy = fit.adequacy, 
        predictivity = predictivities, predictions = predictions)
}
