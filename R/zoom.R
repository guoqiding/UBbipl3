zoom <-
function (x) 
{
    begin <- locator(1)
    usr <- par("usr")
    xrange <- usr[2] - usr[1]
    yrange <- usr[4] - usr[3]
    c(begin$x, begin$x + x * xrange, begin$y, begin$y + x * yrange)
}
