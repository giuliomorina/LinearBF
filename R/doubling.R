doubling_Nacu <- function(n, toss.coin, eps, n0=NULL, max_iter=Inf, verbose=F) {
  #f(p) = min(2p,1-2eps) from Nacu-Peres with Latuszynski envelope method
  doubling_alg(toss.coin, eps, n0, max_iter, verbose)
}