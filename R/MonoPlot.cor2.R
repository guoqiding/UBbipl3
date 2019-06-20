MonoPlot.cor2 <-
function (X = Flotation.data[, -1], as.axes = TRUE, axis.col = "darkgrey", 
    ax.name.col = "black", ax.name.size = 0.65, arrows = TRUE, 
    calib.no = 5, circle = FALSE, dim.plot = 2, exp.factor = 1, 
    n.int = rep(5, ncol(X)), offset = c(0.5, 0.5, 0.5, 0.5), 
    offset.m = rep(0, ncol(X)), pos = "Orthog", pos.m = rep(1, 
        ncol(X)), rotate.degrees = 0, reflect = c(FALSE, "x", 
        "y")) 
{
    n <- nrow(X)
    p <- ncol(X)
    Rmat <- cor(X)
    Rmat2 <- Rmat^2
    centmat <- diag(p) - matrix(1, nrow = p, ncol = p)/p
    PCO.R1 <- svd(0.5 * centmat %*% Rmat %*% centmat)
    PCO.R2 <- svd(0.5 * centmat %*% Rmat2 %*% centmat)
    VSigmaHalf.R <- PCO.R1$v %*% diag(sqrt(PCO.R1$d))
    VSigmaHalf.R2 <- PCO.R2$v %*% diag(sqrt(PCO.R2$d))
    Jmat <- diag(c(rep(1, dim.plot), rep(0, min(n, p) - dim.plot)))
    vars.points.cor1 <- VSigmaHalf.R %*% Jmat
    vars.points.cor2 <- VSigmaHalf.R2 %*% Jmat
    par(pty = "s")
    marker.mat <- vector("list", p)
    plot(vars.points.cor1, xlim = range(vars.points.cor1[, 1] * 
        exp.factor), ylim = range(vars.points.cor1[, 2] * exp.factor), 
        type = "p", pch = 15, col = "red", cex = 1.2, asp = 1, 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:p) text(x = vars.points.cor1[i, 1], y = vars.points.cor1[i, 
        2], label = colnames(Rmat)[i], pos = pos.m[i], offset = offset.m[i], 
        cex = 0.6, xpd = NA)
    points(0, 0, pch = 3)
    if (as.axes) {
        for (i in 1:p) {
            abline(reg = lm(c(0, vars.points.cor1[i, 2]) ~ c(0, 
                vars.points.cor1[i, 1])))
            theta <- atan(coefficients(lm(c(0, vars.points.cor1[i, 
                2]) ~ c(0, vars.points.cor1[i, 1])))[2])
            if (circle) {
                draw.circle.line(r = 1, theta = theta, calib = calib.no, 
                  col = "green", pos.m = pos.m[i], offset.m = offset.m[i])
                draw.circle.line(r = 1, theta = theta + pi, calib = calib.no, 
                  col = "green", pos.m = pos.m[i], offset.m = offset.m[i])
            }
        }
    }
    dev.new()
    par(pty = "s")
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
    rotate.mat <- matrix(c(cos(radns), -sin(radns), sin(radns), 
        cos(radns)), ncol = 2)
    vars.points.cor2 <- vars.points.cor2[, 1:2] %*% rotate.mat %*% 
        reflect.mat
    plot(vars.points.cor2, xlim = range(vars.points.cor2[, 1] * 
        exp.factor), ylim = range(vars.points.cor2[, 2] * exp.factor), 
        type = "p", pch = 15, col = "red", cex = 1.2, asp = 1, 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:p) text(x = vars.points.cor2[i, 1], y = vars.points.cor2[i, 
        2], label = colnames(Rmat)[i], pos = pos.m[i], offset = offset.m[i], 
        cex = 0.6, xpd = NA)
    points(0, 0, pch = 3)
    if (as.axes) {
        for (i in 1:p) {
            abline(reg = lm(c(0, vars.points.cor2[i, 2]) ~ c(0, 
                vars.points.cor2[i, 1])))
            theta <- atan(coefficients(lm(c(0, vars.points.cor2[i, 
                2]) ~ c(0, vars.points.cor2[i, 1])))[2])
            if (circle) {
                draw.circle.line(r = 1, theta = theta, calib = calib.no, 
                  col = "green", pos.m = pos.m[i], offset.m = offset.m[i])
                draw.circle.line(r = 1, theta = theta + pi, calib = calib.no, 
                  col = "green", pos.m = pos.m[i], offset.m = offset.m[i])
            }
        }
    }
    adequacies.R <- diag(((PCO.R1$v[, 1:dim.plot]) %*% t(PCO.R1$v[, 
        1:dim.plot]))/((PCO.R1$v) %*% t(PCO.R1$v)))
    names(adequacies.R) <- colnames(Rmat)
    adequacies.R2 <- diag(((PCO.R2$v[, 1:dim.plot]) %*% t(PCO.R2$v[, 
        1:dim.plot]))/((PCO.R2$v) %*% t(PCO.R2$v)))
    names(adequacies.R2) <- colnames(Rmat2)
    list(cor(X), adequacies.R = adequacies.R, adequacies.R2 = adequacies.R2)
}
