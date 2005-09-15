.onLoad <- function(...) {
   cat("##\n## MSBVAR Package\n")
   cat("## Copyright (C) 2005 Patrick T. Brandt\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0351179, SES-0351205, and SES-0540816)\n##\n")
   require(KernSmooth, quietly=TRUE)
   require(xtable, quietly=TRUE)
}
