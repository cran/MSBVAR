.onLoad <- function(...) {
   cat("##\n## MSBVAR Package v.0.2.2\n")
   cat("## 2007-01-14\n")
   cat("## Copyright (C) 2005-2007 Patrick T. Brandt\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0351179, SES-0351205, and SES-0540816)\n##\n")
   require(KernSmooth, quietly=TRUE)
   require(xtable, quietly=TRUE)
   require(coda, quietly=TRUE)
   require(methods, quietly=TRUE)
}
