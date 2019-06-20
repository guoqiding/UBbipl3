boxplot.NJ <-
function (data = NULL, means = FALSE, notch = FALSE, notchfill = TRUE, 
    confcol = c("grey", "red"), names = NULL, width = NULL, varwidth = FALSE, 
    outline = TRUE, notch.frac = 0.5, log = "", border = par("fg"), 
    col = par("bg"), medcol = "white", pars = list(whisklty = 1), 
    frame.plot = axes, horizontal = FALSE, add = FALSE, at = NULL, 
    show.names = NULL, medld = 2, pch.means = "*", ...) 
{
    if (length(levels(factor(data[, 2]))) == 1) {
        mean.vec <- mean(data[, 1])
        z <- boxplot(data[, 1], plot = FALSE, notch = notch, 
            names = names)
    }
    else {
        mean.vec <- NULL
        for (i in levels(factor(data[, 2]))) {
            mean.vec[i] <- mean(data[data[, 2] == i, ][, 1])
        }
        z <- boxplot(formula = data[, 1] ~ data[, 2], plot = FALSE, 
            notch = notch, names = names)
    }
    pars <- c(pars, list(...))
    if (!missing(col)) 
        warning("argument 'col' is deprecated in favour of 'boxfill'")
    bplt <- function(x, mean.x, wid, stats, out, conf, notch, 
        notchfill, confcol, xlog, i, pch.means, medcol) {
        ok <- TRUE
        if (!any(is.na(stats))) {
            xP <- if (xlog) 
                function(x, w) x * exp(w)
            else function(x, w) x + w
            wid <- wid/2
            if (notch) {
                ok <- stats[2] <= conf[1] && conf[2] <= stats[4]
                xx <- xP(x, wid * c(-1, 1, 1, notch.frac, 1, 
                  1, -1, -1, -notch.frac, -1))
                yy <- c(stats[c(2, 2)], conf[1], stats[3], conf[2], 
                  stats[c(4, 4)], conf[2], stats[3], conf[1])
                if (notchfill) {
                  notch.x <- xP(x, wid * c(-notch.frac, -1, 1, 
                    notch.frac, -notch.frac, -1, 1, notch.frac, 
                    1))
                  notch.y <- c(stats[3], conf[1], conf[1], stats[3], 
                    stats[3], conf[2], conf[2], stats[3])
                  polygon(notch.x[1:4], notch.y[1:4], col = confcol[2])
                  polygon(notch.x[5:8], notch.y[5:8], col = confcol[2])
                  x1 <- xP(x, wid * c(-1, -1, 1, 1, -1, -1, 1, 
                    1))
                  y1 <- c(stats[2], conf[1], conf[1], stats[2], 
                    stats[4], conf[2], conf[2], stats[4])
                  if (notch.y[2] < y1[1]) {
                    polygon(x1[1:4], y1[1:4])
                  }
                  if (notch.y[6] > y1[5]) {
                    polygon(x1[5:8], y1[5:8])
                  }
                  if (notch.y[2] > y1[1]) {
                    polygon(x1[1:4], y1[1:4], col = confcol[1])
                  }
                  if (notch.y[6] < y1[5]) {
                    polygon(x1[5:8], y1[5:8], col = confcol[1])
                  }
                }
            }
            else {
                xx <- xP(x, wid * c(-1, 1, 1, -1))
                yy <- stats[c(2, 2, 4, 4)]
            }
            if (!notch) 
                notch.frac <- 1
            wntch <- notch.frac * wid
            xypolygon(xx, yy, lty = "blank", col = boxfill[i])
            xysegments(xP(x, -wntch), stats[3], xP(x, +wntch), 
                stats[3], lty = medlty[i], lwd = medld, col = medcol[i])
            xypoints(x, stats[3], pch = medpch[i], cex = medcex[i], 
                col = medcol[i], bg = medbg[i])
            xysegments(rep.int(x, 2), stats[c(1, 5)], rep.int(x, 
                2), stats[c(2, 4)], lty = whisklty[i], lwd = whisklwd[i], 
                col = whiskcol[i])
            xysegments(rep.int(xP(x, -wid * staplewex), 2), stats[c(1, 
                5)], rep.int(xP(x, +wid * staplewex), 2), stats[c(1, 
                5)], lty = staplelty[i], lwd = staplelwd[i], 
                col = staplecol[i])
            xypolygon(xx, yy, lty = boxlty[i], lwd = boxlwd[i], 
                border = boxcol[i])
            if ((nout <- length(out)) > 0) {
                xysegments(rep(x - wid * outwex, nout), out, 
                  rep(x + wid * outwex, nout), out, lty = outlty[i], 
                  lwd = outlwd[i], col = outcol[i])
                xypoints(rep.int(x, nout), out, pch = outpch[i], 
                  cex = outcex[i], col = outcol[i], bg = outbg[i])
            }
            if (any(inf <- !is.finite(out))) {
                warning(sprintf(ngettext(length(unique(out[inf])), 
                  "Outlier (%s) in boxplot %d is not drawn", 
                  "Outliers (%s) in boxplot %d are not drawn"), 
                  paste(unique(out[inf]), collapse = ", "), x), 
                  domain = NA)
            }
        }
        if (means) {
            points(x, mean.x, pch = pch.means)
        }
        return(ok)
    }
    if (!is.list(z) || 0 == (n <- length(z$n))) 
        stop("invalid first argument")
    if (is.null(at)) 
        at <- 1:n
    else if (length(at) != n) 
        stop("'at' must have same length as 'z$n', i.e. ", n)
    if (is.null(z$out)) 
        z$out <- numeric()
    if (is.null(z$group) || !outline) 
        z$group <- integer()
    if (is.null(pars$ylim)) 
        ylim <- range(z$stats[is.finite(z$stats)], z$out[is.finite(z$out)], 
            if (notch) z$conf[is.finite(z$conf)])
    else {
        ylim <- pars$ylim
        pars$ylim <- NULL
    }
    if (is.null(pars$xlim)) 
        xlim <- c(0.5, n + 0.5)
    else {
        xlim <- pars$xlim
        pars$xlim <- NULL
    }
    if (length(border) == 0) 
        border <- par("fg")
    if (!add) {
        plot.new()
        if (horizontal) 
            plot.window(ylim = xlim, xlim = ylim, log = log, 
                xaxs = pars$yaxs)
        else plot.window(xlim = xlim, ylim = ylim, log = log, 
            yaxs = pars$yaxs)
    }
    xlog <- (par("ylog") && horizontal) || (par("xlog") && !horizontal)
    pcycle <- function(p, def1, def2 = NULL) rep(if (length(p)) p else if (length(def1)) def1 else def2, 
        length.out = n)
    p <- function(sym) pars[[sym, exact = TRUE]]
    boxlty <- pcycle(pars$boxlty, pars$lty, par("lty"))
    boxlwd <- pcycle(pars$boxlwd, pars$lwd, par("lwd"))
    boxcol <- pcycle(pars$boxcol, border)
    boxfill <- pcycle(pars$boxfill, col)
    boxwex <- pcycle(pars$boxwex, 0.8 * {
        if (n <= 1) 
            1
        else quantile(diff(sort(if (xlog) 
            log(at)
        else at)), 0.1)
    })
    medlty <- pcycle(pars$medlty, p("lty"), par("lty"))
    medlwd <- pcycle(pars$medlwd, 3 * p("lwd"), 3 * par("lwd"))
    medpch <- pcycle(pars$medpch, NA_integer_)
    medcex <- pcycle(pars$medcex, pars$cex, par("cex"))
    medcol <- pcycle(pars$medcol, medcol)
    medbg <- pcycle(pars$medbg, p("bg"), par("bg"))
    whisklty <- pcycle(pars$whisklty, p("lty"), "dashed")
    whisklwd <- pcycle(pars$whisklwd, p("lwd"), par("lwd"))
    whiskcol <- pcycle(pars$whiskcol, border)
    staplelty <- pcycle(pars$staplelty, p("lty"), par("lty"))
    staplelwd <- pcycle(pars$staplelwd, p("lwd"), par("lwd"))
    staplecol <- pcycle(pars$staplecol, border)
    staplewex <- pcycle(pars$staplewex, 0.5)
    outlty <- pcycle(pars$outlty, "blank")
    outlwd <- pcycle(pars$outlwd, p("lwd"), par("lwd"))
    outpch <- pcycle(pars$outpch, p("pch"), par("pch"))
    outcex <- pcycle(pars$outcex, p("cex"), par("cex"))
    outcol <- pcycle(pars$outcol, border)
    outbg <- pcycle(pars$outbg, p("bg"), par("bg"))
    outwex <- pcycle(pars$outwex, 0.5)
    width <- if (!is.null(width)) {
        if (length(width) != n | any(is.na(width)) | any(width <= 
            0)) 
            stop("invalid boxplot widths")
        boxwex * width/max(width)
    }
    else if (varwidth) 
        boxwex * sqrt(z$n/max(z$n))
    else if (n == 1) 
        0.5 * boxwex
    else rep.int(boxwex, n)
    if (horizontal) {
        xypoints <- function(x, y, ...) points(y, x, ...)
        xypolygon <- function(x, y, ...) polygon(y, x, ...)
        xysegments <- function(x0, y0, x1, y1, ...) segments(y0, 
            x0, y1, x1, ...)
    }
    else {
        xypoints <- points
        xypolygon <- polygon
        xysegments <- segments
    }
    ok <- TRUE
    for (i in 1:n) ok <- ok & bplt(at[i], mean.x = mean.vec[i], 
        wid = width[i], stats = z$stats[, i], out = z$out[z$group == 
            i], conf = z$conf[, i], notch = notch, notchfill = notchfill, 
        confcol = confcol, medcol = medcol, pch.means = pch.means, 
        xlog = xlog, i = i)
    if (!ok) 
        warning("Some totches went outside hinges ('box'): maybe set notch=FALSE")
    axes <- is.null(pars$axes)
    if (!axes) {
        axes <- pars$axes
        pars$axes <- NULL
    }
    if (axes) {
        ax.pars <- pars[names(pars) %in% c("xaxt", "yaxt", "las", 
            "cex.axis", "col.axis", "format")]
        if (is.null(show.names)) 
            show.names <- n > 1
        if (show.names) 
            do.call("axis", c(list(side = 1 + horizontal, at = at, 
                labels = z$names), ax.pars))
        do.call("Axis", c(list(x = z$stats, side = 2 - horizontal), 
            ax.pars))
    }
    do.call("title", pars[names(pars) %in% c("main", "cex.main", 
        "col.main", "sub", "cex.sub", "col.sub", "xlab", "ylab", 
        "cex.lab", "col.lab")])
    if (frame.plot) 
        box()
    invisible(at)
}
