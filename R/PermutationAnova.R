PermutationAnova <-
function (X, group.vec, dist = c("Pythagoras", "Clark", "SqrtL1"), 
    n.perm = 2) 
{
    dist <- dist[1]
    if (dist == "Clark") 
        dist.expr <- function(xik, xjk) {
            ((xik - xjk)/(xik + xjk))^2
        }
    if (dist == "Pythagoras") 
        dist.expr <- function(xik, xjk) {
            (xik - xjk)^2
        }
    if (dist == "SqrtL1") 
        dist.expr <- function(xik, xjk) {
            abs(xik - xjk)
        }
    D.fun <- function(mat) {
        n <- nrow(mat)
        D <- matrix(0, nrow = n, ncol = n)
        for (i in 1:(n - 1)) for (j in (i + 1):n) D[i, j] <- -0.5 * 
            sum(dist.expr(mat[i, ], mat[j, ]))
        D + t(D)
    }
    dvec.fun <- function(vec, mat) {
        if (length(vec) != ncol(mat)) 
            stop("Incorrect sized vector")
        n <- nrow(mat)
        outvec <- rep(NA, n)
        for (i in 1:n) outvec[i] <- -0.5 * sum(dist.expr(vec, 
            mat[i, ]))
        outvec
    }
    Gmat <- indmat(group.vec)
    groupnames <- colnames(Gmat)
    Nmat <- t(Gmat) %*% Gmat
    nvec <- diag(Nmat)
    numberOfGroups <- length(nvec)
    Dmat <- D.fun(X)
    AOD.perm <- function(A, N, G, g, n.groups) {
        D <- D.fun(A)
        D.bar <- solve(N) %*% t(G) %*% D %*% G %*% solve(N)
        mat.1 <- matrix(1, nrow = g, ncol = g)
        Deltabar <- D.bar - (diag(diag(D.bar)) %*% mat.1)/2 - 
            (mat.1 %*% diag(diag(D.bar)))/2
        T.ss <- -sum(D)/sum(n.groups)
        B.ss <- -(t(n.groups) %*% Deltabar %*% n.groups)/sum(n.groups)
        return(c(T.ss, B.ss))
    }
    out.null <- AOD.perm(A = X, N = Nmat, G = Gmat, g = numberOfGroups, 
        n.groups = nvec)
    out.perm <- matrix(NA, nrow = n.perm, ncol = 2)
    for (i in (1:n.perm)) {
        vec <- 1:nrow(X)
        num <- nrow(X)
        while (num > nvec[numberOfGroups]) {
            k <- sample(1:num, size = 1)
            temp <- vec[k]
            vec[k] <- vec[num]
            vec[num] <- temp
            num <- num - 1
        }
        out.perm[i, ] <- AOD.perm(A = X[vec, ], N = Nmat, G = Gmat, 
            g = numberOfGroups, n.groups = nvec)
    }
    list(TandB.ss.null = out.null, out.perm.TandB = out.perm, 
        asl = sum(out.perm[, 2] > out.null[2])/n.perm)
}
