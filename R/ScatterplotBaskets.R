ScatterplotBaskets <-
function (x = Headdimensions.data[, 5], y = Headdimensions.data[, 
    6], n = 36, pp = c(0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 
    0.995, 0.998), lty = 1, colr = "blue", lwd = 1) 
{
    if (length(x) != length(y)) 
        stop("x and y are not of the same length! \n")
    expectile <- function(x, a) {
        u <- mean(x)
        for (it in 1:20) {
            w <- a * (x > u) + (1 - a) * (x <= u)
            u0 <- u
            u <- sum(w * x)/sum(w)
            if (u == u0) 
                break
        }
        return(u)
    }
    exp_conts <- function(x, y, p = 0.9, n = 12) {
        a <- sn <- cs <- rep(0, n + 2)
        a0 <- 0
        for (j in 1:(n + 2)) {
            h <- 2 * pi * j/n
            sn[j] <- sin(h)
            cs[j] <- cos(h)
            z <- cs[j] * x + sn[j] * y
            a[j] <- expectile(z, p)
        }
        qx <- qy <- rep(0, n + 1)
        S <- cbind(sn, cs)
        for (j in 1:(n + 1)) {
            r <- j:(j + 1)
            u <- solve(S[r, ], a[r])
            qy[j] <- u[1]
            qx[j] <- u[2]
        }
        return(list(x = qx, y = qy))
    }
    mx <- mean(x)
    sx <- sd(x)
    my <- mean(y)
    sy <- sd(y)
    xs <- (x - mx)/sx
    ys <- (y - my)/sy
    m <- length(x)
    par(mfrow = c(1, 1))
    n <- n
    phis <- (1:(n + 1)) * 2 * pi/n
    sn <- sin(phis)
    cs <- cos(phis)
    uqs <- 0 * phis
    pp <- pp
    np <- length(pp)
    cols <- rainbow(np)
    j <- 0
    for (p in pp) {
        j <- j + 1
        k <- 0
        for (phi in phis) {
            k <- k + 1
            u <- xs * cs[k] + ys * sn[k]
            v <- -xs * sn[k] + ys * cs[k]
            uq <- expectile(u, p)
            xq <- uq * cs[k]
            yq <- uq * sn[k]
            uqs[k] <- uq
            g <- 10 * c(-1, 1)
        }
        qp <- exp_conts(xs, ys, p = p, n = n)
        lines(qp$x * sx + mx, qp$y * sy + my, pch = 15, col = colr, 
            lty = lty, lwd = lwd)
    }
    points(mx, my, pch = 3, cex = 2)
}
