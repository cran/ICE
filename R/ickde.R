"ickde" <-
function (I, h, f, m, n.iterations = 10, x1, xm, right.limit = 10000) 
{
    if (missing(x1)) 
        x1 <- min(I) - 4 * h
    if (missing(xm)) 
        xm <- max(c(I[, 1], I[I[, 2] < right.limit, 2])) + 4 * 
            h
    if (missing(m)) 
        m <- length(f)
    x <- seq(x1, xm, length = m)
    xdiff <- x[2] - x[1]
    if (missing(f)) 
        f <- rep(1, m)/(m * xdiff)
    n <- dim(I)[1]
    left <- I[, 1]
    right <- I[, 2]
    numgrid <- m
    gridpts <- x
    f0 <- f
    niter <- n.iterations
    f1 <- f
    z <- .Fortran("ickde", as.integer(n), as.double(left), as.double(right), 
        as.integer(numgrid), as.double(gridpts), as.double(f0), 
        as.double(h), as.integer(niter), as.double(f1), PACKAGE = "ICE")
    names(z) <- c("n", "left", "right", "numgrid", "x", "f0", 
        "h", "niter", "y")
    z.out <- list(x = z$x, y = z$y)
    class(z.out) <- "IC"
    z.out
}
.First.lib <- function(lib,pkg) {
 library.dynam("ICE",pkg,lib)
}


