doubling_Nacu <-
  function(n,
           toss.coin,
           eps,
           n0 = NULL,
           max_iter = Inf,
           verbose = F) {
    if (is.infinite(max_iter)) {
      max_iter <- -1
    }
    #f(p) = min(2p,1-2eps) from Nacu-Peres with Latuszynski envelope method
    doubling_alg(toss.coin, eps, n0, max_iter, verbose)
  }

list_to_matrix <- function(res_list) {
  res <-
    matrix(unlist(res_list),
           ncol = length(res_list),
           byrow = FALSE)
  rownames(res) <- c("res", "num_tosses")
  return(res)
}

logistic <- function(n, C, p, as_list = F) {
  # f(p) = Cp/(1+Cp)
  # The function is taken from "Designing perfect simulation algorithms
  # using local correctness" by M. Huber (2019), in an iterative form rather than
  # recursive. This is equivalent to the 2-coin algorithm proposed in
  # "Barker's algorithm for Bayesian inference with intractable likelihoods"
  # by F. Goncalves et al. (2017)
  res_list <- lapply(1:n, function(iter) {
    logistic_cpp(C=C, p=p)
  })
  if (as_list) {
    return(res_list)
  } else {
    return(list_to_matrix(res_list))
  }
}

doubling_Huber_2019_rec <-
  function(n,
           C,
           i,
           eps,
           p,
           num_tosses = 0,
           as_list = F) {
    # f(p) = (Cp)^i, for Cp < 1-eps
    # The function is taken from "Designing perfect simulation algorithms
    # using local correctness" by M. Huber (2019).
    # The recursive form, although correct, can sometimes go to deep in the
    # number of recursion, thus giving raise to an error.
    res_list <- lapply(1:n, function(iter) {
      if (i == 0) {
        return(list(res = 1, num_tosses = num_tosses))
      }
      else if (i > 3.55 / eps) {
        beta <- (1 - eps / 2) / (1 - eps)
        B1 <- sample(c(1, 0),
                     size = 1,
                     prob = c(beta ^ (-i), 1 - beta ^ (-i)))
        if (B1 == 0) {
          return(list(res = 0, num_tosses = num_tosses))
        } else {
          return(
            doubling_Huber_2019_rec(
              n = 1,
              C = beta * C,
              i = i,
              eps = eps / 2,
              p = p,
              num_tosses = num_tosses,
              as_list = T
            )
          )
        }
      } else {
        B2 <- logistic(
          n = 1,
          C = C,
          p = p,
          as_list = T
        )[[1]]
        num_tosses <- num_tosses + B2$num_tosses
        return(
          doubling_Huber_2019_rec(
            n = 1,
            C = C,
            i = i + 1 - 2 * B2$res,
            eps = eps,
            p = p,
            num_tosses = num_tosses,
            as_list = T
          )
        )
      }
    })
    if (as_list) {
      return(res_list)
    } else {
      return(list_to_matrix(res_list))
    }
  }

doubling_Huber_2019_iter <-
  function(n,
           C,
           i,
           eps,
           p,
           as_list = F) {
  # f(p) = (Cp)^i, for Cp < 1-eps
  # The function is taken from "Designing perfect simulation algorithms
  # using local correctness" by M. Huber (2019).
  # This is a translation in an iterative form.
  res_list <- lapply(1:n, function(iter) {
    doubling_Huber_2019_iter_cpp(C=C,i=i,eps=eps,p=p)
  })
  if (as_list) {
    return(res_list)
  } else {
    return(list_to_matrix(res_list))
  }
  }

doubling_Huber_2017 <- function() {
  # f(p) = Cp, for Cp < 1-eps
  # The function is taken from "Optimal Linear Bernoulli Factories for 
  # Small Mean Problems" by M. Huber (2017).
}

doubling_Huber_2014 <- function() {
  # f(p) = Cp, for Cp < 1-eps
  # The function is taken from "Nearly optimal Bernoulli factories for 
  # linear functions" by M. Huber (2014).
}

doubling_Thomas_2018_Nacu_Peres <- function() {
  # f(p) = ???
  # The function is taken from "A Practical Implementation of the
  # Bernoulli Factory", which makes use of Latuszynski envelope method
  # and the envelopes provided by Nacu-Peres for twice differentiable functions
}

doubling_Thomas_2018_New_Envelopes <- function() {
  # f(p) = ???
  # The function is taken from "A Practical Implementation of the
  # Bernoulli Factory", which makes use of Latuszynski envelope method
  # and the new envelopes proposed for twice differentiable functions
}