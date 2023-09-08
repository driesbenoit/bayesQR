# bayesQR: Bayesian Quantile Regression

Bayesian quantile regression using the asymmetric Laplace distribution, both continuous as well as binary dependent variables are supported.
The package consists of implementations of the methods of [Yu & Moyeed (2001)](https://doi.org/10.1016/S0167-7152(01)00124-9), [Benoit & Van den Poel (2012)](https://doi.org/10.1002/jae.1216) and [Al-Hamzawi, Yu & Benoit (2012)](https://doi.org/10.1177/1471082X1101200304).
To speed up the calculations, the Markov Chain Monte Carlo core of all algorithms is programmed in Fortran and called from R.

# Cite
To cite bayesQR in publications use:

Benoit DF, Van den Poel D (2017). "bayesQR: A Bayesian Approach to Quantile Regression." Journal of Statistical Software, 76(7), 1–32. (https://doi.org/10.18637/jss.v076.i07).

Corresponding BibTeX entry:
````
@Article{,
  title = {{bayesQR}: A Bayesian Approach to Quantile Regression},
  author = {Dries F. Benoit and Dirk {Van den Poel}},
  journal = {Journal of Statistical Software},
  year = {2017},
  volume = {76},
  number = {7},
  pages = {1--32},
  doi = {10.18637/jss.v076.i07},
}
````

# CRAN
[](https://CRAN.R-project.org/package=bayesQR)
