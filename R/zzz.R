# Start up and unload functions for the MSBVAR package
# See package LICENSE for more details
#

.onAttach <- function(...) {
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
    cat("##\n## MSBVAR Package v.0.5.0\n")
    cat("## 2010-07-01\n")
    cat("## Copyright (C) 2005-", yr, ", Patrick T. Brandt\n", sep="")
    cat("## Written by Patrick T. Brandt\n")
    cat("##\n## Support provided by the U.S. National Science Foundation\n")
    cat("## (Grants SES-0351179, SES-0351205, SES-0540816, and SES-0921051)\n##\n")
    require(KernSmooth, quietly=TRUE)
    require(xtable, quietly=TRUE)
    require(coda, quietly=TRUE)
    require(bit, quietly=TRUE)
    require(mvtnorm, quietly=TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("MSBVAR", libpath)
}
