expand.interval <-
function (interval, left, right) 
{
    if (interval[2] < interval[1]) 
        stop("Interval must be in the form: (minimum, maximum). \n")
    lower <- interval[1]
    upper <- interval[2]
    range <- upper - lower
    if (left == 0) 
        lower <- lower
    else lower <- lower - left * range
    if (right == 0) 
        upper <- upper
    else upper <- upper + right * range
    if (upper < lower) 
        stop("The new interval's upper boundary is smaller than its lower boundary! \n")
    c(lower, upper)
}
