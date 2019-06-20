integral.low <-
function (xvec, yvec, num) 
{
    x <- xvec[1:num]
    y <- yvec[1:num]
    if (num%%2 == 1) {
        h <- (x[num] - x[1])/(num - 1)
        som <- y[1] + y[num] + 2 * sum(y[(2:num)%%2 == 1]) + 
            4 * sum(y[(2:num)%%2 == 0])
        return(area = (som * h)/3)
    }
    else {
        h <- (x[num] - x[1])/(num - 1)
        som <- y[1] + y[num] + 4 * sum(y[(2:num)%%2 == 1]) + 
            2 * sum(y[(2:num)%%2 == 0])
        return(area = (som * h)/3)
    }
}
