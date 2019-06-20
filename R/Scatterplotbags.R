Scatterplotbags <-
function (Z, alpha = 0.5, point.char = 16, point.char.size = 1, 
    point.colour = "black", colours = "red", max.num, Tukey.median = TRUE, 
    line.width = 1, c.hull.n = 10, lty = 1, expand = 0.02) 
{
    if (length(alpha) > 1) {
        if (any(alpha - c(alpha[-1], 0) < 0)) 
            stop("alpha is either a scalar value or a vector with decreasing elements \n")
    }
    if (length(alpha) != length(colours)) 
        stop("colours and alpha must be vectors of the same size")
    .bags.plot <- function(Z, usr, c.hull.n, alpha1, max.num, 
        point.char, point.char.size, line.width, point.colour, 
        colours, lty) {
        x <- Z[, 1]
        y <- Z[, 2]
        if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
            if (!is.numeric(x)) 
                stop(message = "x is not a numeric dataframe or vector.")
        }
        if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
            if ((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x, 
                data.class) == "numeric"))) 
                stop(message = "x is not a numeric dataframe or vector.")
        }
        x <- as.matrix(x)
        if (dim(x)[2] != 1) 
            stop(message = "x is not a vector.")
        if (is.vector(y) || (is.matrix(y) && !is.data.frame(y))) {
            if (!is.numeric(y)) 
                stop(message = "y is not a numeric dataframe or vector.")
        }
        if ((!is.matrix(y) && !is.vector(y)) || is.data.frame(y)) {
            if ((!is.data.frame(y) && !is.numeric(y)) || (!all(sapply(y, 
                data.class) == "numeric"))) 
                stop(message = "y is not a numeric dataframe or vector.")
        }
        y <- as.matrix(y)
        if (dim(y)[2] != 1) 
            stop(message = "y is not a vector.")
        if (nrow(x) != nrow(y)) 
            stop(message = "x and y should have the same length!")
        na.x <- !is.finite(x)
        na.y <- !is.finite(y)
        ok <- !(na.x | na.y)
        x <- x[ok, , drop = FALSE]
        y <- y[ok, , drop = FALSE]
        n <- nrow(x)
        if (length(x) == 0) 
            stop(message = "All observations have missing values")
        if (n == 1) 
            stop(message = "The sample size should be at least two!")
        dimny <- dimnames(y)[[1]]
        if (length(dimny) == 0) 
            dimny <- 1:n
        if (n * alpha1/100 < c.hull.n) {
            hull <- chull(x, y)
            points(x, y, pch = point.char, col = point.colour)
            polygon(x[hull], y[hull], density = 0, col = colours, 
                lwd = line.width, lty = lty)
            stop(paste("Samples too few for ", alpha1/100, "bag.\n", 
                sep = ""))
        }
        else {
            storage.mode(x) <- "double"
            storage.mode(y) <- "double"
            interpx <- rep(0, 2 * n)
            storage.mode(interpx) <- "double"
            interpy <- rep(0, 2 * n)
            storage.mode(interpy) <- "double"
            datatyp <- matrix(0, n, 3)
            storage.mode(datatyp) <- "double"
            datatyp2 <- matrix(0, n, 2)
            storage.mode(datatyp2) <- "double"
            pxpy <- matrix(0, n, 3)
            storage.mode(pxpy) <- "double"
            whisk <- 2
            abagplot.uit <- .Fortran("abagplot", as.integer(n), 
                as.integer(alpha1), x, y, as.integer(whisk), 
                tukm = double(2), interpx = interpx, interpy = interpy, 
                num = as.integer(0), datatyp = datatyp, indoutl = integer(n), 
                datatyp2 = datatyp2, pxpy = pxpy, boxpl = as.integer(0), 
                nointer = as.integer(0), PACKAGE = "UBbipl")
            tukmedian <- abagplot.uit$tukm
            x.vec <- abagplot.uit$interpx
            y.vec <- abagplot.uit$interpy
            if (all(x.vec == 0) & all(y.vec == 0)) 
                stop(message = " x and y both null vectors")
            nie.nul <- !((x.vec == 0) & (y.vec == 0))
            if (Tukey.median) {
                points(x = tukmedian[1], y = tukmedian[2], pch = 16, 
                  cex = 1, col = "red")
            }
            polygon(x.vec[nie.nul], y.vec[nie.nul], density = 0, 
                col = colours, lwd = line.width, lty = as.vector(lty))
        }
    }
    extension.x <- diff(range(Z[, 1])) * expand * c(-1, 1)
    extension.y <- diff(range(Z[, 2])) * expand * c(-1, 1)
    plot(Z[, 1], Z[, 2], xlab = dimnames(Z)[[2]][1], ylab = dimnames(Z)[[2]][2], 
        type = "p", xlim = range(Z[, 1]) + extension.x, ylim = range(Z[, 
            2]) + extension.y, xaxs = "i", yaxs = "i", asp = 1, 
        pch = point.char, cex = point.char.size)
    usr <- par("usr")
    for (i in 1:length(alpha)) {
        alpha.1 <- alpha[i]
        colours.1 <- colours[i]
        if (alpha.1 < 0 | alpha.1 > 0.99) 
            stop(message = "alpha not to be negative or larger than 0.99")
        alpha.entered <- alpha.1
        alpha.1 <- round(alpha.1, digits = 2)
        if (abs(alpha.entered - alpha.1) > 0) 
            cat("alpha has been rounded to ", alpha.1, "\n")
        alpha1 <- 100 * alpha.1
        .bags.plot(Z = Z, usr = usr, alpha1 = alpha1, max.num = max.num, 
            line.width = line.width, point.char = point.char, 
            point.char.size = point.char.size, c.hull.n = c.hull.n, 
            point.colour = point.colour, colours = colours.1, 
            lty = lty)
    }
}
