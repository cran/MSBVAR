.onLoad <- function(...) {
   require(KernSmooth, quietly=TRUE)
   require(xtable, quietly=TRUE)
   require(coda, quietly=TRUE)
   require(bit, quietly=TRUE)
   require(mvtnorm, quietly=TRUE)
   cat("##\n## MSBVAR Package v.0.4.0\n")
   cat("## 2009-06-12\n")
   cat("## Copyright (C) 2005-2009 Patrick T. Brandt\n")
   cat("## Written by Patrick T. Brandt\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0351179, SES-0351205, and SES-0540816)\n##\n")

}

