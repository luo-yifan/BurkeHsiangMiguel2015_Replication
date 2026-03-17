
# cgmreg.R
# R equivalent of cgmreg.ado (Cameron-Gelbach-Miller multi-way clustered SEs)
#
# The fixest package implements CGM multi-way clustered standard errors natively.
# This wrapper provides a familiar interface for code originally using Stata's cgmreg.
#
# Usage:
#   source("script/cgmreg.R")
#   fit <- cgmreg(y ~ x1 + x2 | fe_var, cluster_vars = c("iso_id", "year"), data = df)
#   summary(fit)   # shows multi-way clustered SEs
#   coef(fit)      # coefficients
#   vcov(fit)      # CGM variance-covariance matrix

library(fixest)

#' Cameron-Gelbach-Miller multi-way clustered regression
#'
#' @param formula     fixest-style formula: lhs ~ rhs | fixed_effects
#' @param cluster_vars  character vector of clustering variable names
#' @param data        data frame or data table
#' @param ...         additional arguments passed to feols()
#' @return  feols object with multi-way clustered standard errors
cgmreg <- function(formula, cluster_vars, data, ...) {
  cluster_fml <- as.formula(paste("~", paste(cluster_vars, collapse = " + ")))
  feols(formula, cluster = cluster_fml, data = data, ...)
}
