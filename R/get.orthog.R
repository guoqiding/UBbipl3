get.orthog <-
function (direction, width = 0.05) 
{
    mM <- par3d("modelMatrix")
    pM <- par3d("projMatrix")
    v.mat <- rbind(cbind(0, direction), 1)
    u.mat <- pM %*% mM %*% v.mat
    e.mat <- apply(u.mat, 2, function(x) x[1:2]/x[4])
    orthog.e.mat <- cbind(0, c((e.mat[2, 2] - e.mat[2, 1])/(e.mat[1, 
        2] - e.mat[1, 1]), -1))
    length <- sqrt(sum(orthog.e.mat[, 2]^2))
    orthog.e.mat <- orthog.e.mat * width/length
    orthog.u.mat <- rbind(orthog.e.mat, 0, 1)
    orthog.v.mat <- solve(pM %*% mM) %*% orthog.u.mat
    orthog.line.mat <- apply(orthog.v.mat, 2, function(x) x[1:3]/x[4])
    t(scale(t(orthog.line.mat), center = T, scale = F))
}
