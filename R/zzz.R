.onLoad <- function(...) {
   cat("##\n## MSBVAR Package v.0.1.2\n")
   cat("## Copyright (C) 2005-2006 Patrick T. Brandt\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0351179, SES-0351205, and SES-0540816)\n##\n")
   require(KernSmooth, quietly=TRUE)
   require(xtable, quietly=TRUE)
}
