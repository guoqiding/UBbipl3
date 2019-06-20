Scatterdensity <-
function (Z, expand = 0.02, point.char = 16, point.char.size = 1, 
    point.col = "black", smooth.n = 100, cuts.density = 50, colours.density = c("green", 
        "yellow", "red"), draw.densitycontours = FALSE, layout.heights = c(140, 
        25), line.size = 2.5, between.columns = -1, dens.ax.cex = 0.6, 
    dens.ax.tcl = -0.2, dens.mgp = c(0, -0.25, 0), bandwidth = NULL, 
    parplotmar = rep(4, 4), parlegendmar.dens = c(4, 5, 0, 5)) 
{
    old.par <- par(no.readonly = TRUE)
    graphics.off()
    layout(matrix(1:2, ncol = 1), heights = layout.heights)
    par(pty = "s", mar = parplotmar)
    on.exit(par(old.par))
    require(MASS)
    extension.x <- diff(range(Z[, 1])) * expand * c(-1, 1)
    extension.y <- diff(range(Z[, 2])) * expand * c(-1, 1)
    plot(Z[, 1], Z[, 2], xlab = dimnames(Z)[[2]][1], ylab = dimnames(Z)[[2]][2], 
        type = "n", xlim = range(Z[, 1]) + extension.x, ylim = range(Z[, 
            2]) + extension.y, xaxs = "i", yaxs = "i", asp = 1, 
        pch = point.char, cex = point.char.size)
    usr <- par("usr")
    par(new = TRUE, pty = "s")
    if (!is.null(bandwidth)) 
        ff1 <- kde2d(Z[, 1], Z[, 2], n = smooth.n, lims = usr, 
            h = bandwidth)
    else ff1 <- kde2d(Z[, 1], Z[, 2], n = smooth.n, lims = usr)
    levels <- pretty(range(ff1$z), n = cuts.density)
    col.use <- colorRampPalette(colours.density)
    col.use <- col.use(length(levels) - 1)
    image(ff1, xlim = usr[1:2], ylim = usr[3:4], breaks = levels, 
        xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", asp = 1, col = col.use)
    par(new = TRUE, pty = "s")
    plot(Z[, 1], Z[, 2], xlab = dimnames(Z)[[2]][1], ylab = dimnames(Z)[[2]][2], 
        type = "p", xlim = usr[1:2], ylim = usr[3:4], xaxs = "i", 
        yaxs = "i", asp = 1, pch = point.char, col = point.col, 
        cex = point.char.size)
    par(pty = "m", mar = parlegendmar.dens)
    plot(range(levels), y = 1:2, ylim = c(10, 100), xaxs = "i", 
        yaxs = "i", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
        type = "n", frame.plot = FALSE)
    rect(xleft = levels[-length(levels)], ybottom = 10, xright = levels[-1], 
        ytop = 50, col = col.use, border = FALSE)
    axis(side = 1, at = pretty(levels, n = 8), labels = pretty(levels, 
        n = 8), line = 0, cex.axis = dens.ax.cex, mgp = dens.mgp, 
        tcl = dens.ax.tcl, las = 0)
}
