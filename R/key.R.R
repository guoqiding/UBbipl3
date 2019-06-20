key.R <-
function (x, y, ..., title = "", align = TRUE, background = 0, 
    border = 0, between = 2, corner = c(missing(x), 1), divide = 3, 
    transparent = FALSE, cex = par("cex"), cex.title = 1.5 * 
        max(cex), col = par("col"), lty = par("lty"), lwd = par("lwd"), 
    font = par("font"), pch = par("pch"), adj = 0, type = "l", 
    size = 5, columns = 1, between.columns = 3, angle = 0, density = -1, 
    space = NULL, text.width.multiplier = 1) 
{
    check.types <- function(x, type, types) {
        n <- length(types)
        match.type <- pmatch(type, types, duplicates.ok = TRUE)
        if (any(is.na(match.type))) 
            stop(paste(x, " must be \"", paste(types[-n], collapse = "\", \""), 
                ", or \"", types[n], "\"", sep = ""))
        types[match.type]
    }
    x <- x
    y <- y
    oldxpd <- par("xpd")
    on.exit(par(xpd = oldxpd), TRUE)
    rest <- list(...)
    colnames <- names(rest)
    for (i in (1:length(colnames))[colnames == "text"]) {
        if (!is.list(rest[[i]])) 
            rest[[i]] <- list(rest[[i]])
    }
    actions <- c("points", "lines", "rectangles", "text")
    colnames <- check.types("key components", colnames, actions)
    nrows <- max(sapply(unlist(rest, recursive = FALSE), length))
    ncols <- length(colnames)
    if (!missing(cex)) {
        oldcex <- par("cex")
        on.exit(par(cex = oldcex), TRUE)
        par(cex = mean(cex))
    }
    cx <- par("cxy")[1]
    cy <- par("cxy")[2]
    cx1 <- cx/par("cex")
    cy1 <- cy/par("cex")
    replen <- function(a, b, n) rep(if (is.null(a)) b else a, 
        length = n)
    for (j in seq(ncols)) {
        this <- rest[[j]]
        this$cex <- replen(this$cex, cex, nrows)
        this$size <- replen(this$size, size, nrows)
        this$type <- replen(this$type, type, nrows)
        this$density <- replen(this$density, density, nrows)
        this$angle <- replen(this$angle, angle, nrows)
        this$col <- replen(this$col, col, nrows)
        this$lty <- replen(this$lty, lty, nrows)
        this$lwd <- replen(this$lwd, lwd, nrows)
        this$adj <- replen(this$adj, adj, nrows)
        this$font <- replen(this$font, font, nrows)
        this$pch <- replen(this$pch, pch, nrows)
        rest[[j]] <- this
    }
    text.adj <- width <- height <- matrix(0, nrows, ncols)
    between <- rep(between, length = ncols)
    for (j in seq(ncols)) {
        this <- rest[[j]]
        for (i in seq(nrows)) {
            switch(colnames[j], points = {
                sx <- sy <- this$cex[i]
            }, lines = {
                sx <- this$size[i]
                sy <- this$cex[i]
            }, rectangles = {
                sx <- this$size[i]
                sy <- this$cex[i]
            }, text = {
                sx <- nchar(this[[1]][i]) * this$cex[i] * text.width.multiplier
                sy <- this$cex[i]
                text.adj[i, j] <- this$adj[i]
            })
            width[i, j] <- sx * cx1 + between[j] * cx
            height[i, j] <- sy * cy1
        }
    }
    if (columns != 1) {
        slice <- function(x, p) {
            m <- nrow(x)
            n <- ncol(x)
            if (m%%p != 0) 
                x <- rbind(x, matrix(0, p - m%%p, n))
            q <- nrow(x)/p
            dim(x) <- c(q, p, n)
            x <- aperm(x, c(1, 3, 2))
            dim(x) <- c(q, n * p)
            x
        }
        width[, ncols] <- width[, ncols] + cx * between.columns
        width <- slice(width, columns)
        nc <- ncol(width)
        width[, nc] <- width[, nc] - cx * between.columns
        height <- slice(height, columns)
        text.adj <- slice(text.adj, columns)
    }
    nc <- ncol(width)
    nr <- nrow(width)
    if (align) 
        for (j in seq(nc)) width[, j] <- max(width[, j])
    xpos <- ypos <- matrix(0, nr, nc)
    for (j in seq(length = nc - 1)) xpos[, j + 1] <- xpos[, j] + 
        width[, j]
    xmax <- x + max(xpos + width)
    Between <- rep(between, each = nrow(xpos))
    xpos <- xpos + x + cx * Between * 0.5
    i <- text.adj != 0
    if (any(i)) 
        xpos[i] <- (xpos + text.adj * (width - Between * cx))[i]
    for (i in seq(nr)) height[i, ] <- max(height[i, ])
    ypos[, ] <- y - cumsum(height[, 1])
    if (nchar(title)) 
        ypos <- ypos - cy1 * cex.title
    ymin <- min(ypos)
    title.excess <- x + cex.title * nchar(title) * cx1 * text.width.multiplier - 
        xmax
    if (title.excess > 0) 
        xmax <- xmax + title.excess
    x.offset <- (x - xmax) * corner[1]
    xpos <- xpos + x.offset + max(0, title.excess/2)
    y.offset <- (y - ymin) * (1 - corner[2])
    ypos <- ypos + y.offset + 0.5 * height
    if (background == 0) 
        border <- as.numeric(border)
    if (!transparent) 
        polygon(c(x, xmax, xmax, x) + x.offset, c(y, y, ymin, 
            ymin) + y.offset + 0, col = background, border = border)
    if (nchar(title)) 
        text((x + xmax)/2 + x.offset, y + y.offset - cex.title/2 * 
            cy1, title, adj = 0.5, cex = cex.title)
    if (columns != 1) {
        restack <- function(x, p) {
            n <- ncol(x)/p
            q <- nrow(x)
            dim(x) <- c(q, n, p)
            x <- aperm(x, c(1, 3, 2))
            dim(x) <- c(p * q, n)
            x
        }
        xpos <- restack(xpos, columns)[1:nrows, , drop = FALSE]
        ypos <- restack(ypos, columns)[1:nrows, , drop = FALSE]
    }
    for (j in seq(ncols)) {
        this <- rest[[j]]
        for (i in seq(nrows)) {
            switch(colnames[j], points = {
                points(xpos[i, j], ypos[i, j], cex = this$cex[i], 
                  col = this$col[i], font = this$font[i], pch = this$pch[[i]])
            }, lines = {
                if (this$type[i] != "p") {
                  lines(xpos[i, j] + seq(0, 1, length = divide) * 
                    cx1 * this$size[i], rep(ypos[i, j], divide), 
                    cex = this$cex[i], lwd = this$lwd[i], type = this$type[i], 
                    lty = this$lty[i], pch = this$pch[[i]], font = this$font[i], 
                    col = this$col[i])
                  if (this$type[i] == "b" || this$type[i] == 
                    "o") points(xpos[i, j] + seq(0, 1, length = divide) * 
                    cx1 * this$size[i], rep(ypos[i, j], divide), 
                    cex = this$cex[i], lwd = this$lwd[i], type = "p", 
                    lty = 1, pch = this$pch[[i]], font = this$font[i], 
                    col = this$col[i])
                } else points(xpos[i, j] + 0.5 * cx1 * this$size[i], 
                  ypos[i, j], cex = this$cex[i], lwd = this$lwd[i], 
                  type = this$type[i], lty = this$lty[i], pch = this$pch[[i]], 
                  font = this$font[i], col = this$col[i])
            }, rectangles = {
                polygon(xpos[i, j] + c(0, rep(this$size[i] * 
                  cx1, 2), 0), ypos[i, j] + cy1 * this$cex[i] * 
                  c(-0.5, -0.5, 0.5, 0.5), col = this$col[i], 
                  angle = this$angle[i], density = this$density[i], 
                  border = FALSE)
            }, text = {
                text(xpos[i, j], ypos[i, j], this[[1]][i], adj = this$adj[i], 
                  cex = this$cex[i], col = this$col[i], font = this$font[i])
            })
        }
    }
}
