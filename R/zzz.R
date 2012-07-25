# Start up functions for the MSBVAR package
#
# 2012-01-19 : updated to use new loading mechanism.
# 2012-05-20 : automated build date

msg <- function(...)
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
    cat("##\n## MSBVAR Package v.0.7-1\n")
    cat("## Build date: ", date(), "\n")
    cat("## Copyright (C) 2005-", yr, ", Patrick T. Brandt\n", sep="")
    cat("## Written by Patrick T. Brandt\n")
    cat("##\n## Support provided by the U.S. National Science Foundation\n")
    cat("## (Grants SES-0351179, SES-0351205, SES-0540816, and SES-0921051)\n##\n")
}

.onAttach <- function(...)
{

    msg()
}

.onUnload <- function(libpath) {
    library.dynam.unload("MSBVAR", libpath)
}
