\name{IsraelPalestineConflict}
\docType{data}
\alias{IsraelPalestineConflict}
\title{Weekly Goldstein Scaled Israeli-Palestinian Conflict Data, 1979-2003}
\description{
   This data set gives Goldstein scaled totals that summarize the
   Israeli-Palestinian conflict from April 1979 - December 2003.  These
   are dyadic or directed actions from one party towards the other.  Data
   are weekly starting April 15, 1979.  Positive values indicate
   cooperation, negative values indicate aggression.  These are a subset
   of the Levant dataset from the Kansas / Penn State Event Data System Levant
   dataset.  The data are from Reuters and AFP news sources and encoded
   into the World Event Interaction Survey (WEIS) coding system and
   Goldstein scalings using the Penn State Event Data System TABARI
   program.  Source data can be found on the KEDS site below.
 }
 \usage{data(IsraelPalestineConflict)}
 \format{A matrix containing 1278 observations.  Column one, "i2p", is
   the Israeli actions towards the Palestinians and column two, "p2i" is
   the Palestinian actions towards the Israelis.}
 
 \source{Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
   Time Series Modeling and the Study of Politics: Theory Testing,
   Forecasting, and Policy Analysis"
   \emph{Political Analysis} 14(1):1-36.
 }
 
 \references{
   Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
   Time Series Modeling and the Study of Politics: Theory Testing,
   Forecasting, and Policy Analysis"
   \emph{Political Analysis} 14(1):1-36.

   Goldstein, Joshua. S. 1992. "A Conflict-Cooperation Scale for WEIS
   Event Data" \emph{Journal of Conflict Resolution}. 36:369-385.
   
   Computational Event Data Project \url{http://eventdata.parusanalytics.com}
}
 
\keyword{datasets}
