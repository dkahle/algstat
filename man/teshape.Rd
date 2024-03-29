% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teshape.r
\name{teshape}
\alias{teshape}
\title{Interconvert data structures}
\usage{
teshape(data, out = c("freq", "tab", "raw"), var)
}
\arguments{
\item{data}{a data frame or array.}

\item{out}{the output format, see examples.}

\item{var}{the name of the frequency variable in the dataset, if not
\code{"freq"}.}
}
\value{
a matrix containing the Markov basis as its columns (for easy
addition to tables)
}
\description{
Interconvert an array, a raw data frame, and frequency distribution
data.frame.
}
\details{
Multivariate categorical data can be represented in several ways. Three comon
ways are : a contingency table, a data frame of raw observations (1 row = 1
subject), and a long data frame with a variable containing the counts in the
contingency table.
}
\examples{

# converting a talbe to a data frame
(tab <- structure(
  array(1:8, c(2,2,2)),
  .Dimnames = list(
    A = c("a1", "a2"),
    B = c("b1", "b2"),
    C = c("c1", "c2")
  )
))

teshape(tab, "freq")
teshape(tab, "raw") # nrow = sum(1:8)


# converting a summarized data frame into a table or raw data frame
(data <- teshape(tab, "freq"))
teshape(data, "tab")
stats::xtabs(freq ~ ., data = data)
teshape(data, "tab") == stats::xtabs(freq ~ ., data = data)
teshape(data, "raw")


# converting a raw data frame into a table or summarized data frame
(data <- teshape(tab, "raw"))
teshape(data, "tab")
teshape(data, "freq")


}
\seealso{
\code{\link[stats:xtabs]{stats::xtabs()}}
}
