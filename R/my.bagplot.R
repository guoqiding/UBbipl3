my.bagplot <-
function (x, y, plotinbag = FALSE, plotoutbag = FALSE, ident = FALSE, 
    drawfence = FALSE, drawloop = TRUE, whisk = 1, truncxmin = NULL, 
    truncxmax = NULL, truncymin = NULL, truncymax = NULL, colr.1 = 0, 
    colr.2 = 14, colr.3 = 16, max.num = 2500, ...) 
{
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
    if (n > max.num) {
        n <- max.num
        aa <- sample(1:max.num, replace = FALSE)
        x <- x[aa]
        y <- y[aa]
    }
    dimny <- dimnames(y)[[1]]
    if (length(dimny) == 0) 
        dimny <- 1:n
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    interpol <- matrix(0, n * 2, 2)
    storage.mode(interpol) <- "double"
    datatyp <- matrix(0, n, 3)
    storage.mode(datatyp) <- "double"
    datatyp2 <- matrix(0, n, 2)
    storage.mode(datatyp2) <- "double"
    pxpy <- matrix(0, n, 3)
    storage.mode(pxpy) <- "double"
    res <- .Fortran("bagplot", as.integer(n), x, y, as.integer(whisk), 
        tukm = double(2), interpol = interpol, num = as.integer(0), 
        datatyp = datatyp, indoutl = integer(n), datatyp2 = datatyp2, 
        pxpy = pxpy, boxpl = as.integer(0), nointer = as.integer(0), 
        PACKAGE = "BookBipl")
    nucleus <- 0
    if (exists("windows")) {
        hulpcol <- matrix(1, 2)
        hulpcol[1] <- colr.1
        hulpcol[2] <- colr.2
    }
    if (res$boxpl != 0) {
        name <- par()$pty
        par(pty = "s")
        eps <- 1e-06
        z <- cbind(x, y)
        dimnames(z) <- list(dimny, NULL)
        h <- res$interpol[1:2, 1:2]
        if (abs(h[2, 1] - h[1, 1]) > eps) {
            rc <- atan((h[2, 2] - h[1, 2])/(h[2, 1] - h[1, 1]))
            const <- h[1, 2] - (h[1, 1] * tan(rc))
            vert <- 0
        }
        else vert <- 1
        if (vert == 1) {
            outl <- z[abs(z[, 1] - h[1, 1]) > eps, , drop = FALSE]
            z <- z[abs(z[, 1] - h[1, 1]) <= eps, , drop = FALSE]
        }
        else {
            outl <- z[abs(z[, 2] - h[1, 2] - tan(rc) * (z[, 1] - 
                h[1, 1])) > eps, , drop = FALSE]
            z <- z[abs(z[, 2] - h[1, 2] - tan(rc) * (z[, 1] - 
                h[1, 1])) <= eps, , drop = FALSE]
        }
        if (vert == 1) {
            print("At least 50% of the data lie on the vertical line:", 
                quote = FALSE)
            print(paste(" x=", signif(h[1, 1], digits = 6)), 
                quote = FALSE)
        }
        else {
            if (const >= 0) {
                print("At least 50% of the data lie on the line:", 
                  quote = FALSE)
                print(paste(" y=", signif(tan(rc), digits = 8), 
                  "x+", signif(const, digits = 8)), quote = FALSE)
            }
            else {
                print("At least 50% of the data lie on the line:")
                print(paste(" y=", signif(tan(rc), digits = 8), 
                  "x", signif(const, digits = 8)), quote = FALSE)
            }
        }
        print("Therefore, the bagplot reduces to a univariate boxplot.", 
            quote = FALSE)
        if (vert == 1) {
            hulp <- z[, 1]
            z[, 1] <- z[, 2]
            z[, 2] <- hulp
        }
        else for (i in 1:nrow(z)) {
            hulp <- z[i, 1]
            z[i, 1] <- cos(rc) * hulp + sin(rc) * z[i, 2]
            z[i, 2] <- -sin(rc) * hulp + cos(rc) * z[i, 2]
        }
        q <- quantile(z[, 1], c(0.25, 0.5, 0.75))
        zhelp <- sort(z[, 1])
        rank1 <- floor(n/4)
        q[1] <- zhelp[rank1]
        rank2 <- ceiling((3 * n)/4)
        q[3] <- zhelp[rank2]
        m <- q[2]
        iqr <- q[3] - q[1]
        db <- 0.07 * iqr
        f1 <- q[1] - 1.5 * iqr
        f1 <- m + 3 * (q[1] - m)
        f2 <- q[3] + 1.5 * iqr
        f2 <- m + 3 * (q[3] - m)
        pl <- matrix(0, 4, 2)
        pl[1, 1] <- f1
        pl[1, 2] <- db
        pl[2, 1] <- f1
        pl[2, 2] <- -db
        pl[3, 1] <- f2
        pl[3, 2] <- db
        pl[4, 1] <- f2
        pl[4, 2] <- -db
        if (length(outl) != 0) 
            if (sum((z[, 1] > f2) + (z[, 1] < f1)) != 0) {
                out <- z[(z[, 1] > f2) | (z[, 1] < f1), , drop = FALSE]
                if (vert == 1) {
                  hulp <- out[, 1]
                  out[, 1] <- out[, 2]
                  out[, 2] <- hulp
                }
                else for (i in 1:nrow(out)) {
                  hulp <- out[i, 1]
                  out[i, 1] <- cos(rc) * hulp - sin(rc) * out[i, 
                    2]
                  out[i, 2] <- sin(rc) * hulp + cos(rc) * out[i, 
                    2]
                }
                outl <- rbind(outl, out)
            }
            else if (sum((z[, 1] > f2) + (z[, 1] < f1)) != 0) 
                outl <- z[(z[, 1] > f2) | (z[, 1] < f1), , drop = FALSE]
        if (sum((z[, 1] > f2) + (z[, 1] < f1)) != 0) 
            z <- z[(z[, 1] <= f2) & (z[, 1] >= f1), , drop = FALSE]
        f1 <- min(z[, 1])
        f2 <- max(z[, 1])
        pl <- rbind(pl, matrix(0, 14, 2))
        pl[5, 1] <- f1
        pl[6, 1] <- q[1]
        pl[7, 1] <- q[3]
        pl[8, 1] <- f2
        pl[9, 1] <- q[1]
        pl[9, 2] <- db
        pl[10, 1] <- q[1]
        pl[10, 2] <- -db
        pl[11, 1] <- q[2]
        pl[11, 2] <- db
        pl[12, 1] <- q[2]
        pl[12, 2] <- -db
        pl[13, 1] <- q[3]
        pl[13, 2] <- db
        pl[14, 1] <- q[3]
        pl[14, 2] <- -db
        pl[15, 1] <- q[1]
        pl[15, 2] <- db
        pl[16, 1] <- q[3]
        pl[16, 2] <- db
        pl[17, 1] <- q[1]
        pl[17, 2] <- -db
        pl[18, 1] <- q[3]
        pl[18, 2] <- -db
        if (vert == 1) 
            pl[, 2] <- pl[, 2] + h[1, 1]
        else pl[, 2] <- pl[, 2] + (-sin(rc) * h[1, 1] + cos(rc) * 
            h[1, 2])
        if (vert == 1) {
            hulp <- z[, 1]
            z[, 1] <- z[, 2]
            z[, 2] <- hulp
            hulp <- pl[, 1]
            pl[, 1] <- pl[, 2]
            pl[, 2] <- hulp
        }
        else {
            for (i in 1:nrow(z)) {
                hulp <- z[i, 1]
                z[i, 1] <- cos(rc) * hulp - sin(rc) * z[i, 2]
                z[i, 2] <- sin(rc) * hulp + cos(rc) * z[i, 2]
            }
            for (i in 1:nrow(pl)) {
                hulp <- pl[i, 1]
                pl[i, 1] <- cos(rc) * hulp - sin(rc) * pl[i, 
                  2]
                pl[i, 2] <- sin(rc) * hulp + cos(rc) * pl[i, 
                  2]
            }
        }
        if (length(outl) != 0) {
            pl <- rbind(pl, outl)
        }
        if (length(outl) != 0) {
            points(outl[, 1], outl[, 2], pch = 17, cex = 1)
        }
        for (i in 1:7) lines(pl[(i * 2 + 3):(i * 2 + 4), 1], 
            pl[(i * 2 + 3):(i * 2 + 4), 2], ...)
        title(" ")
        if (ident) {
            plak <- rbind(z, outl)
            lab <- dimnames(plak)[[1]]
            identify(plak[, 1], plak[, 2], lab)
        }
        par(pty = name)
    }
    else {
        print(paste("The coordinates of the Tukey median are (", 
            signif(res$tukm[1], digits = 6), ",", signif(res$tukm[2], 
                digits = 6), ")."), quote = FALSE)
        if (res$nointer == 1) {
            points(x, y, ...)
            points(res$tukm[1], res$tukm[2], pch = 3)
            for (i in 1:n) lines(c(res$tukm[1], x[i]), c(res$tukm[2], 
                y[i]))
            print("Sorry, no bagplot can be calculated, too many points coincide.", 
                quote = FALSE)
        }
        else {
            if (res$nointer == 2) 
                print("There are no points between the bag and the fence, \n                     so no whiskers can be calculated", 
                  quote = FALSE)
            if ((n >= 15) && (drawfence == TRUE)) {
                fence <- res$interpol[1:res$num, ]
                fence[, 1] <- (3 * (res$interpol[1:res$num, 1] - 
                  res$tukm[1])) + res$tukm[1]
                fence[, 2] <- (3 * (res$interpol[1:res$num, 2] - 
                  res$tukm[2])) + res$tukm[2]
                if (!is.null(truncxmin)) 
                  fence[fence[, 1] < truncxmin, 1] <- truncxmin
                if (!is.null(truncxmax)) 
                  fence[fence[, 1] > truncxmax, 1] <- truncxmax
                if (!is.null(truncymin)) 
                  fence[fence[, 2] < truncymin, 2] <- truncymin
                if (!is.null(truncymax)) 
                  fence[fence[, 2] > truncymax, 2] <- truncymax
            }
            if (n < 15) 
                print("The bag is only plotted when there are at least 15 observations.")
            drawdata <- res$datatyp[res$datatyp[, 3] < 3, , drop = FALSE]
            if (plotinbag == FALSE) 
                drawdata <- drawdata[drawdata[, 3] > 1, , drop = FALSE]
            if (plotoutbag == FALSE) 
                drawdata <- res$datatyp[res$datatyp[, 3] < 2, 
                  , drop = FALSE]
            if ((plotinbag == FALSE) && (plotoutbag == FALSE)) 
                drawdata <- matrix(res$tukm, ncol = 2)
            if (length(drawdata) == 0) 
                drawdata <- matrix(res$tukm, ncol = 2)
            loop <- res$datatyp[, 1:2]
            if (n >= 15) {
                if (drawfence == TRUE) 
                  loop <- fence
                outl <- res$datatyp[, 3] == 3
                if (sum(outl) != 0) {
                  outlier <- res$datatyp[outl, , drop = FALSE]
                  plakx <- c(loop[, 1], outlier[, 1])
                  plaky <- c(loop[, 2], outlier[, 2])
                }
                else {
                  plakx <- loop[, 1]
                  plaky <- loop[, 2]
                }
                xlim1 <- min(plakx)
                xlim2 <- max(plakx)
                ylim1 <- min(plaky)
                ylim2 <- max(plaky)
            }
            else {
                drawdata <- res$datatyp[res$datatyp[, 3] < 4, 
                  , drop = FALSE]
                if (plotinbag == FALSE) 
                  drawdata <- drawdata[drawdata[, 3] > 1, , drop = FALSE]
                if (plotoutbag == FALSE) 
                  drawdata <- res$datatyp[res$datatyp[, 3] < 
                    2, , drop = FALSE]
                if ((plotinbag == FALSE) && (plotoutbag == FALSE)) 
                  drawdata <- matrix(res$tukm, ncol = 2)
                if (length(drawdata) == 0) 
                  drawdata <- matrix(res$tukm, ncol = 2)
                xlim1 <- min(res$datatyp[, 1])
                xlim2 <- max(res$datatyp[, 1])
                ylim1 <- min(res$datatyp[, 2])
                ylim2 <- max(res$datatyp[, 2])
            }
            if (n >= 15) {
                if (nucleus) {
                  delta <- 0.1
                  polygon(delta * (res$interpol[1:res$num, 1] - 
                    res$tukm[1]) + res$tukm[1], delta * (res$interpol[1:res$num, 
                    2] - res$tukm[2]) + res$tukm[2], density = -1)
                }
                else points(res$tukm[1], res$tukm[2], pch = 3, 
                  cex = 1)
            }
            if (n >= 15) {
                if (drawloop == TRUE) {
                  chulp <- res$datatyp[res$datatyp[, 3] < 3, 
                    , drop = FALSE]
                  chu <- rbind(chulp[, 1:2], res$interpol[1:res$num, 
                    1:2])
                  ch <- chull(chu[, 1], chu[, 2])
                  print(ch)
                  polygon(chu[ch, 1], chu[ch, 2], density = -1, 
                    lty = 3, col = hulpcol[1], border = FALSE)
                  polygon(res$interpol[1:res$num, 1], res$interpol[1:res$num, 
                    2], density = -1, col = hulpcol[2], border = TRUE)
                  points(res$tukm[1], res$tukm[2], pch = 16, 
                    col = 0, cex = 1)
                  points(res$tukm[1], res$tukm[2], pch = 3)
                }
                if (drawfence == TRUE) 
                  polygon(fence[, 1], fence[, 2], density = 0, 
                    lty = 3)
                if (sum(outl) != 0) {
                  points(outlier[, 1], outlier[, 2], pch = 8, 
                    cex = 1)
                }
            }
            if (n < 15) {
                for (i in 1:n) lines(c(res$tukm[1], res$datatyp[i, 
                  1]), c(res$tukm[2], res$datatyp[i, 2]))
            }
        }
        if (ident) {
            lab <- dimny[res$indoutl]
            identify(res$datatyp[, 1], res$datatyp[, 2], lab)
        }
    }
    invisible()
}
