colour.fun <-
function (R, G, B) 
{
    hexvec <- c(1:9, "A", "B", "C", "D", "E", "F")
    red1 <- floor(R/16)
    red2 <- R - floor(R/16) * 16
    if (red1 == 0) 
        R1 <- "0"
    else R1 <- hexvec[red1]
    if (red2 == 0) 
        R2 <- "0"
    else R2 <- hexvec[red2]
    green1 <- floor(G/16)
    green2 <- G - floor(G/16) * 16
    if (green1 == 0) 
        G1 <- "0"
    else G1 <- hexvec[green1]
    if (green2 == 0) 
        G2 <- "0"
    else G2 <- hexvec[green2]
    blue1 <- floor(B/16)
    blue2 <- B - floor(B/16) * 16
    if (blue1 == 0) 
        B1 <- "0"
    else B1 <- hexvec[blue1]
    if (blue2 == 0) 
        B2 <- "0"
    else B2 <- hexvec[blue2]
    paste("#", R1, R2, G1, G2, B1, B2, sep = "")
}
