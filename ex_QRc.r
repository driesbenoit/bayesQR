# Install bayesQR version
install.packages("bayesQR.tar.gz",repos=NULL)

# load lib
library(bayesQR)

# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- cbind(1,X)
y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6*X[,2])
#y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6)

# Initiate plot
## Plot datapoints
plot(X[,2], y, main="", cex=.6, xlab="X")

# Write loop to analyze 5 quantiles
cnt <- 0
for (i in c(.05,.25,.5,.75,.95)) {
    cnt <- cnt+1

    ## Estimate parameters
    out = bayesQR(y~0+X,quantile=i,ndraw=5000)

    ## Add quantile regression lines to the plot (exclude first 500 burn-in draws)
    abline(a=mean(out[[1]]$betadraw[500:5000,1]),b=mean(out[[1]]$betadraw[500:5000,2]),lty=cnt,col=cnt)

}

# Estimate and plot OLS model
outOLS = lm(y~0+X)
abline(outOLS$coef,lty=1,lwd=2,col=6)

# Add legend to plot
legend(x=0,y=max(y),legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")
