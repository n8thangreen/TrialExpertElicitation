% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_funcs.R
\name{prior_theta}
\alias{prior_theta}
\title{Identify parameters of a prior normal distribution for theta~N(mu, sigma2)}
\usage{
prior_theta(pi1, gamma, a, b, margin)
}
\arguments{
\item{pi1}{P(pE > pC) as identified by Day 1 Elicitation Q3}

\item{gamma}{P(pE - pC > -margin) as identified by 1 - Q4}

\item{a, b}{parameters of prior distribution for pC}

\item{margin}{non-inferiority margin for the trial}
}
\value{
vector with param\link{1} = mu (prior mean of theta), param\link{2} = sigma2 (prior variance of theta).
}
\description{
Identify parameters of a prior normal distribution for theta~N(mu, sigma2)
}
