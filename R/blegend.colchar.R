blegend.colchar <-
function (classes, pch.means, pch.samples, pch.samples.size = c(1, 
    1), colours, colours.means, legend.type = c(means = TRUE, 
    samples = TRUE, bags = FALSE), line.type, line.width = 1, 
    line.size = 2.5, columns = 1, between = c(4, -3, 4, -3, 4), 
    between.columns = 1, parlegendmar = c(1, 1, 1, 1), text.width.mult = 1, 
    quality.print = FALSE, QualityOfDisplay = "") 
{
    if (all(legend.type == FALSE)) 
        return(cat("Change legend.type to obtain a legend\n"))
    par(pty = "m", mar = parlegendmar)
    plot(x = c(0, 10), y = c(0, 10), type = "n", axes = FALSE, 
        xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    usr <- par("usr")
    x <- usr[1]
    y <- usr[4]
    if (quality.print) {
        mtext(text = QualityOfDisplay, adj = 0, cex = 0.8, at = 0, 
            side = 3, line = 0.5, las = 0)
    }
    if (legend.type[3]) {
        if (legend.type[1] & legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.means, col = colours.means, 
                  cex = pch.samples.size[1]), points = list(pch = pch.samples, 
                  col = colours, cex = pch.samples.size[1]), 
                lines = list(lty = line.type, col = colours.means, 
                  lwd = line.width), text = list(classes, cex = pch.samples.size[2], 
                  col = colours), border = T, columns = columns, 
                between.columns = between.columns, size = line.size, 
                text.width.multiplier = text.width.mult)
        if (legend.type[1] & !legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.means, col = colours.means, 
                  cex = pch.samples.size[1]), lines = list(lty = line.type, 
                  col = colours.means, lwd = line.width), text = list(classes, 
                  cex = pch.samples.size[2], col = colours), 
                border = T, columns = columns, between.columns = between.columns, 
                size = line.size, text.width.multiplier = text.width.mult)
        if (!legend.type[1] & legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.samples, col = colours, 
                  cex = pch.samples.size[1]), lines = list(lty = line.type, 
                  col = colours.means, lwd = line.width), text = list(classes, 
                  cex = pch.samples.size[2], col = colours), 
                border = T, columns = columns, between.columns = between.columns, 
                size = line.size, text.width.multiplier = text.width.mult)
        if (!legend.type[1] & !legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                lines = list(lty = line.type, col = colours.means, 
                  lwd = line.width), text = list(classes, cex = pch.samples.size[2], 
                  col = colours), border = T, columns = columns, 
                between.columns = between.columns, size = line.size, 
                text.width.multiplier = text.width.mult)
    }
    else {
        if (legend.type[1] & legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.means, col = colours.means, 
                  cex = pch.samples.size[1]), points = list(pch = pch.samples, 
                  col = colours, cex = pch.samples.size[1]), 
                text = list(classes, cex = pch.samples.size[2], 
                  col = colours), border = T, columns = columns, 
                between.columns = between.columns, size = line.size, 
                text.width.multiplier = text.width.mult)
        if (legend.type[1] & !legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.means, col = colours.means, 
                  cex = pch.samples.size[1]), text = list(classes, 
                  cex = pch.samples.size[2], col = colours), 
                border = T, columns = columns, between.columns = between.columns, 
                size = line.size, text.width.multiplier = text.width.mult)
        if (!legend.type[1] & legend.type[2]) 
            key.R(x = x, y = y, corner = c(0, 1), between = between, 
                points = list(pch = pch.samples, col = colours, 
                  cex = pch.samples.size[1]), text = list(classes, 
                  cex = pch.samples.size[2], col = colours), 
                border = T, columns = columns, between.columns = between.columns, 
                size = line.size, text.width.multiplier = text.width.mult)
    }
}
