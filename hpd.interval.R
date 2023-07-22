## vec: vector of MCMC samples of one parameter
hpd.interval <- function(vec, prob = 0.95)
{
    vals <- sort(vec)
    nsamp <- length(vals)
    npar <- 1
    gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
    init <- 1:(nsamp - gap)
    inds <- which.min(vals[init + gap] - vals[init])
    ans <- c(vals[inds], vals[inds + gap])
    attr(ans, "Probability") <- gap/nsamp
    ans
}