const_a_n_Flegal_Morina <- function(n, k, eps, C) {
  M <- C*C*sqrt(2)/(eps*sqrt(exp(1)))
  exp(Re(cgamma(complex(real=n, imaginary=-sqrt(M)), log=T) + 
           cgamma(complex(real=n, imaginary=sqrt(M)), log=T)) -
        2*lgamma(n))
}