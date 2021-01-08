This repository contains code for the content of Chapter 2 of the PhD thesis "Extending the Bernoulli Factory to a Dice Enterprise".

## Linear Bernoulli Factory

Given a coin that lands head with probability <img src="https://render.githubusercontent.com/render/math?math=p">, a linear Bernoulli Factory is an algorithm that produces tosses of a coin that has probability <img src="https://render.githubusercontent.com/render/math?math=Mp"> of landing heads, under the assumption <img src="https://render.githubusercontent.com/render/math?math=Mp \leq 1-\epsilon, \epsilon > 0">.

Several linear Bernoulli Factories are implemented, in particular:

- `doubling_Nacu_Peres` implements via `Rcpp` the doubling algorithm proposed by Nacu and Peres (2005).
- `doubling_Flegal_Herbei` implements the linear Bernoulli Factory proposed by Flegal and Herbei (2012).
- `doubling_Huber_2016` implements the linear Bernoulli Factory proposed by Huber (2016).
- `doubling_Huber_2017` and `small_doubling_Huber_2017` implement the linear Bernoulli Factory algorithms proposed by Huber (2017). The latter further assumes $M < 1/2$.
- `doubling_Huber_2019_rec` implements the linear Bernoulli Factory proposed by Huber (2019) in its original recursive formulation. The function `doubling_Huber_2019_iter` implements the same algorithm in an iterative form so to avoid stack overflow. 

## Comparison of Huber's algorithms

The vignette includes a comparison of the three algorithms proposed by Huber. Results are presented in Section 2.2 of the thesis.

## References

Nacu, S ̧. and Peres, Y. (2005). Fast simulation of new coins from old. Ann. Appl. Probab. 15(1A), 93–115.

Flegal, J. M. and Herbei, R. (2012). Exact sampling for intractable probability distributions via a Bernoulli factory.Electronic Journal of Statistics 6, 10–37.

Huber, M. (2016). Nearly optimal Bernoulli factories for linear functions.Combinatorics, Probability and Computing 25(4), 577–591.

Huber, M. (2017). Optimal Linear Bernoulli Factories for Small Mean Problems.Methodology and Computing in Applied Probability 19(2), 631–645.

Huber, M. (2019). Designing perfect simulation algorithms using local correctness. arXiv preprint arXiv:1907.06748.
