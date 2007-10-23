.onLoad <- function(...) {
   cat("##\n## MSBVAR Package v.0.3.1\n")
   cat("## 2007-10-23\n")
   cat("## Copyright (C) 2005-2007 Patrick T. Brandt\n")
   cat("## Written by Patrick T. Brandt and Justin Appleby\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0351179, SES-0351205, and SES-0540816)\n##\n")
   require(KernSmooth, quietly=TRUE)
   require(xtable, quietly=TRUE)
   require(coda, quietly=TRUE)
}

