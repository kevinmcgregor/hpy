# hpy
Hierarchical Pitman-Yor Gibbs sampler and diversity estimation

# Instructions for installing hpy

Note this package imports gStirling.  Please install before use:
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github("kevinmcgregor/gStirling")

# To install hpy run:
install_github("kevinmcgregor/hpy")

