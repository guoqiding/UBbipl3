tied.names <-
function (cats = c("a", "b", "c", "d", "e", "f", "g"), vals = c(0.1, 
    0.2, 0.2, 0.3, 0.3, 0.3, 0.4)) 
{
    unique.points <- table(vals)
    k <- length(unique.points)
    cat.labs <- rep(NA, k)
    tel <- 1
    for (i in 1:k) {
        if (unique.points[i] == 1) {
            cat.labs[i] <- cats[tel]
            tel <- tel + 1
        }
        else {
            cat.lab.temp <- cats[tel]
            for (j in 1:(unique.points[i] - 1)) {
                cat.lab.temp <- paste(cat.lab.temp, cats[tel + 
                  j], sep = "*")
            }
            cat.labs[i] <- cat.lab.temp
            tel <- tel + cumsum(unique.points[i])
        }
    }
    cat.labs
}
