draw.halfway.lines <-
function (x, y, reverse = FALSE, colour = 1:10, width = rep(2, 
    10), w.factor = 1, ...) 
{
    regr <- coefficients(lm(y ~ x))
    xmin <- par("usr")[1] - 2
    xmax <- par("usr")[2] + 2
    ymin <- regr[1] + regr[2] * xmin
    ymax <- regr[1] + regr[2] * xmax
    x <- c(xmin, x, xmax)
    y <- c(ymin, y, ymax)
    x.half <- (x[-1] + x[-length(x)])/2
    y.half <- (y[-1] + y[-length(y)])/2
    k <- length(x.half)
    for (i in 1:(k - 1)) {
        tel1 <- i
        tel2 <- i + 1
        if (identical(reverse, FALSE)) 
            lines(x = c(x.half[tel1], x.half[tel2]), y = c(y.half[tel1], 
                y.half[tel2]), col = colour[i], lwd = width[i] * 
                w.factor, ...)
        if (identical(reverse, TRUE)) 
            lines(x = c(x.half[tel1], x.half[tel2]), y = c(y.half[tel1], 
                y.half[tel2]), col = colour[k - i], lwd = width[k - 
                i] * w.factor, ...)
    }
}
