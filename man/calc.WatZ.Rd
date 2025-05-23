% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.WatZ.R
\name{calc.WatZ}
\alias{calc.WatZ}
\title{Convert Tanner crab size to weight using the current (2015) relationships}
\usage{
calc.WatZ(
  z,
  sex,
  maturity,
  male = list(immature = list(a = 0.00027, b = 3.022134), mature = list(a = 0.00027, b =
    3.022134), undetermined = list(a = 0.00027, b = 3.022134)),
  female = list(immature = list(a = 0.000562, b = 2.816928), mature = list(a = 0.000441,
    b = 2.898686))
)
}
\arguments{
\item{z}{: vector of sizes (mm CW)}

\item{sex}{: vector of sexes ('MALE' or 'FEMALE')}

\item{maturity}{: vector of maturity states ('IMMATURE', 'MATURE', or 'UNDETERMINED')}

\item{male}{: list of regression coefficients for immature, mature crab (see details)}

\item{female}{: list of regression coefficients for immature, mature crab (see details)}
}
\value{
a vector of weights in g corresponding to the input sizes, sexes, and maturity states
}
\description{
Power-law functions for Tanner crab weight-at-size (g) by sex and maturity
based on 2015 NMFS survey.
}
\details{
'male' and 'female' are lists with named elements 'a' and 'b', whose values
reflect the parameters in the eq:  \eqn{w = a \cdot z^b}

For males, the default values are a=0.00027, b=3.022134 for immature, mature, and undetermined.

For females, the default values are \cr
a=0.000562, b=2.816928 for immature crab \cr
a=0.000441, b=2.898686 for mature crab \cr
For females with undetermined maturity, females > 90 mm CW are considered mature.
}
