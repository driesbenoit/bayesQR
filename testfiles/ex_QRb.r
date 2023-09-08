# Install bayesQR version
install.packages("bayesQR.tar.gz",repos=NULL)

# load lib
library(bayesQR)

# Simulate data from heteroskedastic regression
set.seed(66)
n = 200
X = as.matrix(runif(n=n, min=0, max=10))
ystar = 1 - 1.5*X + rnorm(n=n, mean=0, sd=2*X)
ystar = 1 - 1.5*X + rnorm(n=n, mean=0, sd=1)
y <- as.numeric(ystar>0)

# Estimate parameters
#out <- bayesQR(y~0+X,quantile=.75,ndraw=5000)
out <- bayesQR(y~X,quantile=.5,ndraw=50000,keep=10)

plot(out,plottype="trace")
summary(out)


##############
set.seed(66)
n = 600
X = cbind(1,runif(n=n, min=-5, max=5))
ystar = X%*%c(1,1.5) + rnorm(n=n, mean=0, sd=1)
y <- as.numeric(ystar>0)

out <- bayesQR(y~0+X,quantile=.5,ndraw=50000,keep=10)
plot(out,plottype="trace")
summary(out)
