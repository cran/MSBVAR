\name{BCFdata}
\docType{data}
\alias{BCFdata}
\alias{Y}
\alias{z2}
\title{Subset of Data from Brandt, Colaresi, and Freeman (2008)}
\description{
  This data set in two matrices about the Israeli-Palestinian conflict.
  The first matrix is a set of endogenous variables that gives 1)
  monthly Goldstein scaled means that summarize the
  Israeli-Palestinian conflict from April 1996 - March 2005, 2) average
  Jewish peace index score [0 = no support  to 100=full support] that
  measure Jewish public support of the
  peace process based on polls of Jewish respondents from the Tami
  Steinmetz Center for Peace Research. The conflict measures 
  are dyadic or directed actions from one party towards the other.
  Positive values indicate an average of more cooperation and less
  conflict and negative values indicate an average with more conflict
  than cooperation.  These are a subset
  of the Levant dataset from the Kansas / Penn State / Computational
  Event Data Project Levant dataset.  The data are from AFP news sources
  and encoded into the World Event Interaction Survey (WEIS) coding
  system and Goldstein scalings using the Event Data Project TABARI
  program.  Source data can be found on the site below.
 }
 \usage{data(BCFdata)}
 \format{Two matrices containing 108 observations.  The first matrix "Y"
   is a multiple \code{ts} object of the endogenous series that measure
   the average conflict-cooperation level and the public opinion data.
   This matrix has three columns. Column one, "I2P", is average
   Goldstein scaled Israeli actions towards the Palestinians; column
   two, "P2I" is average Goldstein scaled Palestinian actions towards
   the Israelis; column three is the average Jewish peace index value
   for the month, "JPI".

   The second matrix, "z2"  is a set of nine control variables for
   shifts in the conflict, the Israeli prime ministerial regime, and
   election trends.  The columns of this matrix are 1) a dummy variable
   for the period from the start of the second Intifada to the start of
   the Battle of Jenin (October 2000--April 2002, end of the second
   Intifada). 2) a dummy variable for the post-Battle of Jenin period
   (May 2002--March 2005), 3-5) dummy variables for the
   identities of the Israeli prime ministers in each month (one each
   for Netanyahu, Barak, and Sharon, with Rabin/Peres treated as the
   reference category). 6-9) a separate time counter that starts at the
   value 1 in the month after each Israeli election and increases until
   the time of the next constitutionally mandated election.
} 
\source{Brandt, Patrick T., Michael P. Colaresi and John R. Freeman.
  2008. ``The Dynamics of Reciprocity, Accountability and
  Credibility.'' \emph{Journal of Conflict Resolution}. 52(3): 343-374.
  
  Replication materials at
  \url{http://jcr.sagepub.com/content/52/3/343}
  
 }
 
 \references{

   Goldstein, Joshua. S. 1992. "A Conflict-Cooperation Scale for WEIS
   Event Data" \emph{Journal of Conflict Resolution}. 36:369-385.
   
   Computational Event Data Project \url{http://eventdata.parusanalytics.com/}
}

\keyword{datasets}
