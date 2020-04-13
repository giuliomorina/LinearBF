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

doubling_Huber_2017 <- function(n, C, eps, p, as_list=F) {
  # f(p) = Cp, for Cp < 1-eps
  # The function is taken from "Optimal Linear Bernoulli Factories for 
  # Small Mean Problems" by M. Huber (2017).
  
  # Define auxiliary functions
  A_func <- function(m,C,p) {
    s <- 1
    num_tosses <- 0
    while(s > 0 && s < m) {
      aux <- logistic(n=1, C=C, p=p, as_list = T)[[1]]
      num_tosses <- num_tosses+aux$num_tosses
      B <- aux$res
      s <- s-2*B+1
    }
    if(s == 0) {
      res <- 1
    } else {
      res <- 0
    }
    return(list(res=res, num_tosses=num_tosses))
  }
  High_Power_Logistic_BF <- function(m, beta, C, p) {
    s <- 1
    num_tosses <- 0
    while(s > 0 && s <= m) {
      aux <- logistic(n=1, C=beta*C, p=p, as_list = T)[[1]]
      num_tosses <- num_tosses+aux$num_tosses
      B <- aux$res
      s <- s-2*B-1
    }
    if(s == m+1) {
      res <- 1
    } else {
      res <- 0
    }
    return(list(res=res, num_tosses=num_tosses))
  }
  B_func <- function(eps, m, beta, C, p) {
    res <- -1
    num_tosses <- 0
    while(res == -1) {
      aux <- doubling_Huber_2017(n=1, C=beta*C, eps=1-(1-eps)*beta, p=p, as_list=T)[[1]]
      B1 <- aux$res
      num_tosses <- num_tosses + aux$num_tosses
      if(B1 == 0) {
        res <- 0
      } else {
        aux2 <- High_Power_Logistic_BF(m=m-2, beta=beta, C=C, p=p)
        B2 <- aux2$res
        num_tosses <- num_tosses + aux2$num_tosses
        if(B2 == 1) {
          res <- 1
        } else {
          m <- m-1
        }
      }
    }
    return(list(res=res, num_tosses=num_tosses))
  }
  # Main function
  res_list <- lapply(1:n, function(iter) {
    m <- ceiling(4.5/eps)+1
    beta <- 1+1/(m-1)
    B1 <- A_func(m, beta*C, p)
    num_tosses <- B1$num_tosses
    if(B1$res == 1) {
      B2 <- sample(c(1,0), size=1, prob=c(1/beta,1-1/beta))
      if(B2 == 1) {
        res <- 1
      } else {
        aux <- B_func(eps, m, beta, C, p)
        res <- aux$res
        num_tosses <- num_tosses + aux$num_tosses
      }
    } else {
      res <- 0
    }
    return(list(res=res, num_tosses = num_tosses))
  })
  if (as_list) {
    return(res_list)
  } else {
    return(list_to_matrix(res_list))
  }
}

small_doubling_Huber_2017 <- function(n, C, M, p, as_list=F, doubling_bf=doubling_Huber_2017) {
  # f(p) = Cp, for Cp <= M < 1/2, with known constant M
  # The function is taken from "Optimal Linear Bernoulli Factories for 
  # Small Mean Problems" by M. Huber (2017).
  res_list <- lapply(1:n, function(iter) {
    num_tosses <- 0
    beta <- 1/(1-2*M)
    aux <- logistic(n=1,C=beta*C,p=p, as_list = T)[[1]]
    num_tosses <- num_tosses + aux$num_tosses
    Y <- aux$res
    B <- sample(c(1,0), size=1, prob=c(1/beta, 1-1/beta))
    if(Y != 1) {
      res <- 0
    } else if(Y==1 && B==1) {
      res <- 1
    } else {
      aux2 <- doubling_bf(n=1,C=beta*C/(beta-1), eps=1-M, p=p, as_list=T)[[1]]
      res <- aux2$res
      num_tosses <- aux2$num_tosses
    }
    return(list(res=res, num_tosses=num_tosses))
  })
  if (as_list) {
    return(res_list)
  } else {
    return(list_to_matrix(res_list))
  }
}

doubling_Huber_2014 <- function(n, C, eps, p, as_list=F) {
  # f(p) = Cp, for Cp < 1-eps
  # The function is taken from "Nearly optimal Bernoulli factories for 
  # linear functions" by M. Huber (2014).
  res_list <- lapply(1:n, function(iter) {
    num_tosses <- 0
    gamma <- 0.5
    k <- 2.3/(gamma*eps)
    i <- 1
    eps <- min(eps, 0.644)
    R <- 1
    while(i > 0 && R==1) {
      while(i > 0 && i < k) {
        B <- sample(c(1,0), size = 1, prob = c(p,1-p))
        num_tosses <- num_tosses + 1
        G <- rgeom(n = 1, prob = (C-1)/C) + 1 #Starting from 1
        i <- i-1+(1-B)*G
      }
      if(i >= k) {
        R <- sample(c(1,0), size = 1, prob=c((1+gamma*eps)^(-i), 1-(1+gamma*eps)^(-i)))
        C <- C*(1+gamma*eps)
        eps <- (1-gamma)*eps
        k <- k/(1-gamma)
      }
    }
    if(i == 0) {
      res <- 1
    } else {
      res <- 0
    }
    return(list(res=res, num_tosses=num_tosses))
  })
  if (as_list) {
    return(res_list)
  } else {
    return(list_to_matrix(res_list))
  }
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