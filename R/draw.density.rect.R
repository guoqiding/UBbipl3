draw.density.rect <-
function (levels.rect, col.use, dens.ax.cex, dens.mgp, dens.ax.tcl) 
{
    plot(range(levels.rect), y = 1:2, ylim = c(10, 100), xaxs = "i", 
        yaxs = "i", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
        type = "n", frame.plot = FALSE)
    rect(xleft = levels.rect[-length(levels.rect)], ybottom = 10, 
        xright = levels.rect[-1], ytop = 50, col = col.use, border = FALSE)
    axis(side = 1, at = pretty(levels.rect, n = 8), labels = pretty(levels.rect, 
        n = 8), line = 0, cex.axis = dens.ax.cex, mgp = dens.mgp, 
        tcl = dens.ax.tcl, las = 0)
}
