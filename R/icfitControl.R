`icfitControl`<-function (epsilon = 1e-06, maxit = 10000, initfitOpts=NULL) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, initfitOpts= initfitOpts)
}
