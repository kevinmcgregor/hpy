#' Title
#'
#' @param Y
#' @param tab
#' @param n.s.tab
#' @param pop
#' @param c.top
#' @param d.top
#' @param c.loc
#' @param d.loc
#'
#' @return
#' @export
#'
#' @examples
getSimpsonConditional <- function(Y, tab, n.s.tab, pop, c.top, d.top, c.loc, d.loc) {

  # Current table frequencies
  Y.c <- Y[pop,] # Species frequencies for this populations
  c.tab <- tab[[pop]] # Tables for population "pop"
  n.tab <- NROW(c.tab) # Number of tables in population "pop"
  n.j <- sum(c.tab[,2]) # Number of individuals in population "pop"
  n.tab.sp <- colSums(n.s.tab) # Number of tables of each species across all populations
  tot.tab <- sum(n.tab.sp) # Number of tables overall
  K <- length(n.tab.sp) # Total number of species across all populations

  d1 <- c.loc + n.j
  d2 <- d1 + 1
  d3 <- c.top + tot.tab
  d4 <- d3 + 1

  ca <- (c.loc+n.tab*d.loc)*(c.top+K*d.top)
  no.sum1 <- ca*(1-d.loc)
  no.sum2 <- ca*(c.loc+(n.tab+1)*d.loc)*(1-d.top)

  sp.vec <- c.tab[, 1]
  sum.t1 <- sum((c.tab[,2]-d.loc)*(Y[pop,sp.vec]+1-(n.s.tab[pop,sp.vec]+1)*d.loc))
  sum.t2 <- sum((c.tab[,2]-d.loc)*(n.tab.sp[sp.vec]-d.top))

  sum.k1 <- sum((n.tab.sp-d.top)*(Y[pop,]+1-(n.s.tab[pop,]+1)*d.loc))
  sum.k2 <- sum((n.tab.sp-d.top)*(c.loc+(n.tab+1)*d.loc)*(n.tab.sp+1-d.top))

  d <- sum.t1/(d1*d2) + (c.loc+n.tab*d.loc)*(sum.t2/(d1*d2*d3) + sum.k1/(d1*d2*d3) + sum.k2/(d1*d2*d3*d4)) +
    no.sum1/(d1*d2*d3) + no.sum2/(d1*d2*d3*d4)

  return(1-d)
}

#' Title
#'
#' @param Y Species table: rows are local populations, columns are species
#' @param tab
#' @param n.s.tab
#' @param pop1
#' @param pop2
#' @param c.top
#' @param d.top
#' @param c.loc1
#' @param d.loc1
#' @param c.loc2
#' @param d.loc2
#'
#' @return
#' @export
#'
#' @examples
getSimpsonBetaConditional <- function(Y, tab, n.s.tab, pop1, pop2, c.top, d.top, c.loc1, d.loc1, c.loc2, d.loc2) {
  Y.1 <- Y[pop1,] # Species frequencies for this population
  Y.2 <- Y[pop2,]
  c.tab1 <- tab[[pop1]]
  c.tab2 <- tab[[pop2]]
  n.tab1 <- NROW(c.tab1)
  n.tab2 <- NROW(c.tab2)
  n.1 <- sum(c.tab1[,2])
  n.2 <- sum(c.tab2[,2])
  n.tab.sp <- colSums(n.s.tab) # Number of tables of each species over all populations
  tot.tab <- sum(n.tab.sp) # Number of tables overall
  K <- NCOL(n.tab.sp)

  d1 <- c.loc1 + n.1
  d2 <- c.loc2 + n.2
  d3 <- c.top + tot.tab
  d4 <- d3 + 1

  no.sum <- (c.top+K*d.top)*(c.loc2+n.tab2*d.loc2)*(1-d.top)

  sp.vec <- c.tab1[, 1]
  cd <- (c.tab1[,2]-d.loc1)
  t.sum1 <- sum(cd*(Y.2[sp.vec]-n.s.tab[pop2,sp.vec]*d.loc2))
  t.sum2 <- sum(cd*(n.tab.sp[sp.vec]-d.top))

  nk <- (n.tab.sp-d.top)
  k.sum1 <- sum(nk*(Y.2-n.s.tab[pop2,]*d.loc2))
  k.sum2 <- (c.loc2+n.tab2*d.loc2)*sum(nk*(n.tab.sp+1-d.top))

  d <- t.sum1/(d1*d2) + (c.loc2+n.tab2*d.loc2)*t.sum2/(d1*d2*d3) + (c.loc1+n.tab1*d.loc1)*(
    (no.sum+k.sum2)/(d1*d2*d3*d4) + k.sum1/(d1*d2*d3))
  return(1-d)
}

#' Calculate unconditional Simpson's index for a particular population
#'
#' @param c.top Top-level concentration parameter
#' @param d.top Top-level discount parameter
#' @param c.loc Local-level concentration parameter for the population of interest
#' @param d.loc Local-level discount parameter for the population of interest
#'
#' @return Unconditional Simpson's index
#' @export
#'
#' @examples
getSimpsonUnconditional <- function(c.top, d.top, c.loc, d.loc) {
  1 - (1-d.loc)/(c.loc+1) - (c.loc+d.loc)/(c.loc+1)*(1-d.top)/(c.top+1)
}

#' Calculate unconditional Simpson's index for beta diversity
#'
#' @param c.top Top-level concentration parameter
#' @param d.top Top-level discount parameter
#'
#' @return Unconditional Simpson's index for beta diversity
#' @export
#'
#' @examples
getSimpsonBetaUnconditional <- function(c.top, d.top) {
  1 - (1-d.top)/(c.top+1)
}

