CircularNonLinear.predictions <-
function (ToProject, g = c(0, 0), Biplot.axis, Biplot.markers) 
{
    if (length(ToProject) == 2) 
        ToProject <- as.numeric(ToProject)
    else stop("ToProject must be a 2-element vector \n")
    p.in <- length(Biplot.axis)
    eps <- 1e-08
    Biplot.points.WhereClosestOnAxis.temp <- NULL
    Biplot.points.WhichClosestOnAxis.temp <- NULL
    lambda.vec <- NULL
    lambda.temp <- NULL
    c.vec <- NULL
    d.vec <- NULL
    A.vec <- (ToProject + g[1:2])/2
    r <- sqrt(sum((ToProject - A.vec)^2))
    lambda.func <- function(tempA, tempB) {
        c.vec <<- Biplot.axis[[tempA]][tempB, ]
        d.vec <<- Biplot.axis[[tempA]][tempB + 1, ]
        a.scalar <- c.vec %*% c.vec - 2 * c.vec %*% d.vec + d.vec %*% 
            d.vec
        b.scalar <- 2 * c.vec %*% d.vec - 2 * A.vec %*% c.vec + 
            2 * d.vec %*% A.vec - 2 * d.vec %*% d.vec
        c.scalar <- -2 * d.vec %*% A.vec + d.vec %*% d.vec + 
            A.vec %*% A.vec - r^2
        if (abs(a.scalar) > eps) {
            if (b.scalar^2 - 4 * a.scalar * c.scalar >= -eps) {
                lambda1 <- (-b.scalar + sqrt(max(0, b.scalar^2 - 
                  4 * a.scalar * c.scalar)))/(2 * a.scalar)
                lambda2 <- (-b.scalar - sqrt(max(0, b.scalar^2 - 
                  4 * a.scalar * c.scalar)))/(2 * a.scalar)
            }
            else lambda1 <- lambda2 <- NA
        }
        else lambda1 <- lambda2 <- -c.scalar/b.scalar
        c(lambda1, lambda2)
    }
    interpolate.func <- function(lambda, tempA, tempB) {
        X.vec.temp <- lambda * c.vec + (1 - lambda) * d.vec
        if (temp5 > sqrt(sum((X.vec.temp - ToProject)^2))) {
            X.vec <<- X.vec.temp
            X.vec.marker.label <<- lambda * Biplot.markers[[tempA]][tempB] + 
                (1 - lambda) * Biplot.markers[[tempA]][tempB + 
                  1]
            lambda.temp <<- tempB
            temp5 <<- sqrt(sum((X.vec - ToProject)^2))
        }
    }
    for (tempA in 1:p.in) {
        temp0 <- as.matrix(dist(rbind(A.vec, Biplot.axis[[tempA]])))[-1, 
            1]
        temp1 <- sign(temp0 - r)
        temp2 <- sign(diff(temp1))
        temp3 <- abs(temp2)
        temp4 <- which(temp3 == 1)
        temp5 <- Inf
        X.vec <- NULL
        X.vec.marker.label <- NULL
        temp6 <- which(temp1[-length(temp1)] == 1 & temp3 == 
            0)
        u.mat <- Biplot.axis[[tempA]][temp6 + 1, ] - Biplot.axis[[tempA]][temp6, 
            ]
        v.mat <- sweep(-Biplot.axis[[tempA]][temp6, ], 2, A.vec, 
            "+")
        temp7 <- temp6[temp7oftemp6 <- (rowSums(u.mat * v.mat) >= 
            -eps * rowSums(u.mat^2) & rowSums(u.mat * v.mat) <= 
            (1 + eps) * rowSums(u.mat^2))]
        for (tempB in sort(c(temp4, temp7))) {
            lambda <- lambda.func(tempA, tempB)
            if (!is.na(lambda[1]) && lambda[1] >= -eps && lambda[1] <= 
                1 + eps) 
                interpolate.func(lambda[1], tempA, tempB)
            if (!is.na(lambda[1]) && lambda[2] >= -eps && lambda[2] <= 
                1 + eps) 
                interpolate.func(lambda[2], tempA, tempB)
        }
        if (is.null(X.vec)) {
            X.vec <- c(NA, NA)
            X.vec.marker.label <- NA
        }
        Biplot.points.WhereClosestOnAxis.temp <- matrix(rbind(Biplot.points.WhereClosestOnAxis.temp, 
            X.vec), ncol = 2)
        Biplot.points.WhichClosestOnAxis.temp <- c(Biplot.points.WhichClosestOnAxis.temp, 
            X.vec.marker.label)
        lambda.vec <- c(lambda.vec, lambda.temp)
    }
    list(Markers = Biplot.points.WhichClosestOnAxis.temp, Coordinates = Biplot.points.WhereClosestOnAxis.temp)
}
