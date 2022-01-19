#' Hierarchical Pitman-Yor Gibbs sampler
#'
#' @param Y Matrix of taxon counts (rows are populations, columns are species)
#' @param n.iter Number of MCMC iterations (after burn-in)
#' @param n.burn Number of burn-in iterations
#' @param p.shape Shape parameter for gamma prior on top-level concentration parameter
#' @param p.scale Scale parameter for gamma prior on top-level concentration parameter
#' @param p.shape.pop Vector with length equal to number of populations containing shape parameters for gamma priors on population-level concentration parameter
#' @param p.scale.pop Vector with length equal to number of populations containing scale parameters for gamma priors on population-level concentration parameter
#' @param type Type of hierarchical process: "py" for Pitman-Yor process, "dp" for Dirichlet process
#' @param quiet Set to TRUE to suppress console output
#'
#' @return Returns a list containing the following elements:
#' \itemize{
#'   \item gamma, alpha, theta, sigma - Vectors containing posterior samples for each of the HPY parameters.
#'   \item tab - List with table indicators for each population.
#'   \item t.c - Table counts (and species corresponding to each table)
#'   \item n.tab - Number of tables in each population
#'   \item n.s.tab - Number of tables for a given species in each population
#' }
#' @export
#'
#' @import gStirling
#'
#' @examples
hpySampler <- function(Y, n.iter, n.burn, p.shape, p.scale, p.shape.pop=NULL, p.scale.pop=NULL, type=c("py", "dp"), quiet=FALSE) {
  if (!is.numeric(Y) | !is.matrix(Y) |
      any(Y<0) | any(Y!=floor(Y))) stop("Y must be a numeric matrix of positive counts")

  type <- match.arg(type)

  J <- NROW(Y)
  K <- NCOL(Y) # Number of distinct species in joint sample

  if (is.null(p.shape.pop)) {
    p.shape.pop <- rep(p.shape, J)
  }
  if (is.null(p.scale.pop)) {
    p.scale.pop <- rep(p.scale, J)
  }

  # Initializing PY parameters
  gamma <- 5
  alpha <- ifelse(type=="dp", 0, 0.5)
  theta <- rep(1, J)
  if(type=="dp") {
    sigma <- rep(0, J)
  } else {
    sigma <- rep(0.5, J)
  }

  # Initializing table info
  sp.vec <- vector("list", length=J)
  n <- rowSums(Y) # The number of individuals in each population
  tab <- vector("list", length=J) # List to hold table indicators for each population
  t.c <- vector("list", length=J) # Table counts (and species corresponding to each table)
  n.tab <- rep(0, J) # Number of tables in each population
  #n.s.tab <- matrix(1, J, K) # Number of tables for a given species in each population
  n.s.tab <- 1*(Y>0) # Number of tables for a given species in each population
  for (j in 1:J) {
    # Initially one table for each species
    sp.vec[[j]] <- rep(1:K, Y[j,])
    tab[[j]] <- rep(1:sum(Y[j,]>0), c(Y[j,])[Y[j,]>0])
    t.c[[j]] <- cbind((1:K)[Y[j,]>0], c(Y[j,])[Y[j,]>0])
    n.tab[j] <- NROW(t.c[[j]])
  }

  gamma.s <- rep(0, n.iter)
  alpha.s <- rep(0, n.iter)
  theta.s <- matrix(0, n.iter, J)
  sigma.s <- matrix(0, n.iter, J)
  n.tab.s <- matrix(0, n.iter, J)
  n.s.tab.s <- array(0, dim=c(n.iter, J, K))
  t.c.s <- vector("list", length=n.iter)

  # Loop over MCMC iterations
  for (i in 1:(n.burn+n.iter)) {
    idx <- i - n.burn

    if (!quiet & i==1) cat("Beginning burn-in:", "\n")
    if (!quiet & i==n.burn+1) cat("Beginning sampling:", "\n")
    if (!quiet & i%%100==0) cat(" ", i, "\n")

    # Loop over populations
    for (j in 1:J) {
      #Loop over individuals in a population
      for (p in 1:n[j]){
        sp.cur <- sp.vec[[j]][p]
        # Remove current individual from its table
        t.cur <- tab[[j]][p]
        t.c[[j]][t.cur,2] <- t.c[[j]][t.cur,2] - 1
        if (t.c[[j]][t.cur,2]==0) {
          # Drop table and make adjustments
          t.c[[j]] <- t.c[[j]][-t.cur,]
          tab[[j]][tab[[j]]>t.cur] <- tab[[j]][tab[[j]]>t.cur] - 1
          n.tab[j] <- n.tab[j] - 1
          n.s.tab[j, sp.cur] <- n.s.tab[j, sp.cur] - 1
        }

        # Reassign individual to new or existing table
        npart <- (gamma+sum(n.tab))*(Y[j,sp.cur]-1-n.s.tab[j, sp.cur]*sigma[j])
        dpart <- (theta[j]+n.tab[j]*sigma[j])*(sum(n.s.tab[, sp.cur])-alpha)
        prob.new <- ifelse(n.s.tab[j, sp.cur]==0, 1, 1/(1+npart/dpart))
        is.new <- rbinom(1, 1, prob.new)
        if (is.new) {
          # Allocate new table
          t.c[[j]] <- rbind(t.c[[j]], c(sp.cur, 1))
          n.tab[j] <- n.tab[j] + 1
          n.s.tab[j, sp.cur] <- n.s.tab[j, sp.cur] + 1
          tab[[j]][p] <- n.tab[j]
        } else {
          wh.c.sp <- which(t.c[[j]][,1]==sp.cur)
          wh.t.new <- sample(1:length(wh.c.sp), 1, prob=(t.c[[j]][wh.c.sp,2]-sigma[j]))
          wh.t <- wh.c.sp[wh.t.new]
          t.c[[j]][wh.t, 2] <- t.c[[j]][wh.t, 2] + 1
          tab[[j]][p] <- wh.t
        }
      }


      if (type=="py") {
        sigma[j] <- gStirling::sampDsct(sigma[j], K, n.tab[j],Y[j,,drop=FALSE], n.s.tab[j,,drop=FALSE], theta[j])
      }
      theta[j] <- gStirling::sampConc(theta[j], n[j], n.tab[j], p.shape.pop[j], p.scale.pop[j], sigma[j])
    }

    for (j in 1:J) {
      if (i>n.burn) {
        theta.s[idx,j] <- theta[j]
        sigma.s[idx,j] <- sigma[j]
        n.tab.s[idx,j] <- n.tab[j]
      }
    }

    # Sample top-level PY parameters
    # Discount
    sp.tab.totals <- matrix(apply(n.s.tab, 2, sum), ncol=K)
    if (type=="py") {
      alpha <- gStirling::sampDsct(alpha, K, K, sp.tab.totals, t(rep(1, K)), gamma) # Chris
    }
    # Concentration
    gamma <- gStirling::sampConc(gamma, sum(n.tab), K, p.shape, p.scale, alpha)

    if (i>n.burn) {
      gamma.s[idx] <- gamma
      alpha.s[idx] <- alpha
      n.s.tab.s[idx,,] <- n.s.tab
      t.c.s[[idx]] <- t.c
    }
  }

  return(list(gamma=gamma.s, alpha=alpha.s, theta=theta.s, sigma=sigma.s,
              tab=tab, t.c=t.c.s, n.tab=n.tab.s, n.s.tab=n.s.tab.s))
}
