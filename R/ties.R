ties <-
function (d, delta, method = c("primary", "secondary")) 
{
    method <- method[1]
    if (!(identical(method, "primary") | identical(method, "secondary"))) 
        stop("method must be either 'primary' or 'secondary' ")
    if (is.null(names(delta))) 
        labels <- 1:length(delta)
    else labels <- names(delta)
    pairs <- data.frame(d, delta, labels)
    if (identical(method, "primary")) 
        pairs <- pairs[order(d, delta), ]
    if (identical(method, "secondary")) 
        pairs <- pairs[order(d), ]
    z <- as.vector(pairs[, 2])
    names(z) <- pairs[, 3]
    while (!all(order(z) == (1:length(z)))) {
        block.vec <- rep(1, length(z))
        block.num <- 1
        for (i in 2:length(z)) {
            if (z[i] > z[i - 1]) 
                block.num <- block.num + 1
            block.vec[i] <- block.num
        }
        for (i in 1:max(block.vec)) z[block.vec == i] <- mean(z[block.vec == 
            i])
    }
    z
}
