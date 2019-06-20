GabrielPCABiplot <-
function (datmat = scale(Ocotea.data[, -(1:2)]), sleep = 0.05) 
{
    anim <- function(sleep = 0.05, upoints, vpoints, xlims = c(-10, 
        10), ylims = c(-10, 10)) {
        alphavals <- c(seq(from = 0, to = 1, length = 50), seq(from = 1, 
            to = 0, length = 50))
        for (i in alphavals) {
            plot.new()
            par(pty = "s")
            plot.window(xlims, ylims, asp = 1)
            axis(side = 1, at = pretty(xlims), labels = rep("", 
                length(pretty(xlims))), line = 1)
            axis(side = 2, at = pretty(ylims), labels = rep("", 
                length(pretty(ylims))), line = 1)
            points(0, 0, pch = 3, cex = 2, col = "black")
            points(upoints %*% (diag(uit$d[1:2]^i)), pch = 15, 
                col = "green")
            text(upoints %*% (diag(uit$d[1:2]^i)), labels = dimnames(datmat)[[1]], 
                cex = 0.6)
            points(vpoints %*% (diag(uit$d[1:2]^(1 - i))), pch = 16, 
                col = "red")
            text(vpoints %*% (diag(uit$d[1:2]^(1 - i))), labels = dimnames(datmat)[[2]], 
                cex = 0.6)
            for (j in 1:nrow(vpoints)) {
                lines(c(0, (vpoints %*% (diag(uit$d[1:2]^(1 - 
                  i))))[j, 1]), c(0, (vpoints %*% (diag(uit$d[1:2]^(1 - 
                  i))))[j, 2]))
                lines(c(0, -(vpoints %*% (diag(uit$d[1:2]^(1 - 
                  i))))[j, 1]), c(0, -(vpoints %*% (diag(uit$d[1:2]^(1 - 
                  i))))[j, 2]))
            }
            par(usr = c(0, 1, 0, 1))
            text(x = c(0.1, 0.9), y = c(0.05, 0.05), labels = c("alpha = 0.0", 
                "alpha = 1.0"), cex = 0.8)
            lines(x = c(0.15, 0.15 + i * 0.78), y = c(0.01, 0.01), 
                lwd = 6, col = "blue")
            Sys.sleep(sleep)
        }
        Recall(sleep, upoints, vpoints, xlims, ylims)
    }
    uit <- svd(datmat)
    u2 <- uit$u[, 1:2]
    v2 <- uit$v[, 1:2]
    u2alpha <- diff(range(u2 %*% (diag(uit$d[1:2]))))
    v2alpha <- diff(range(v2 %*% (diag(uit$d[1:2]))))
    if (u2alpha > v2alpha) 
        lims <- range(u2 %*% (diag(uit$d[1:2])))
    else lims <- range(v2 %*% (diag(uit$d[1:2])))
    print(uit)
    cat(paste("Quality = ", sum(uit$d[1:2]^2)/sum(uit$d^2)), 
        "\n")
    anim(sleep = sleep, upoints = u2, vpoints = v2, xlims = lims, 
        ylims = lims)
}
