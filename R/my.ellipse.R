my.ellipse <-
function (mu, sigma.vec, cormat, minvec, colour = colour.fun(153, 
    153, 255), level = 0.9, alpha = 0.5) 
{
    mu <- mu - minvec
    mu <- mu[c(1, 3, 2)]
    sigma.vec <- sigma.vec[c(1, 3, 2)]
    cormat <- cormat[, c(1, 3, 2)]
    cormat <- cormat[c(1, 3, 2), ]
    plot3d(ellipse3d(cormat, scale = sigma.vec, centre = mu, 
        level = level), add = TRUE, col = colour, alpha = alpha)
}
