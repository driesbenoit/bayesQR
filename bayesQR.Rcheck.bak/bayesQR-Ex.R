pkgname <- "bayesQR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bayesQR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bayesQR")
### * bayesQR

flush(stderr()); flush(stdout())

### Name: bayesQR
### Title: Bayesian quantile regression
### Aliases: bayesQR

### ** Examples

# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Estimate series of quantile regressions with adaptive lasso
out <- bayesQR(y~X, quantile=c(.05,.25,.5,.75,.95), alasso=TRUE, ndraw=5000)

# Initiate plot
## Plot datapoints
plot(X, y, main="", cex=.6, xlab="X")

## Add quantile regression lines to the plot (exclude first 500 burn-in draws)
sum <- summary(out, burnin=500)
for (i in 1:length(sum)){
  abline(a=sum[[i]]$betadraw[1,1],b=sum[[i]]$betadraw[2,1],lty=i,col=i)
}

# Estimate and plot OLS model
outOLS = lm(y~X)
abline(outOLS,lty=1,lwd=2,col=6)

# Add legend to plot
legend(x=0,y=max(y),legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
       lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")



cleanEx()
nameEx("plot.bayesQR")
### * plot.bayesQR

flush(stderr()); flush(stdout())

### Name: plot.bayesQR
### Title: Produce quantile plots or traceplots with 'plot.bayesQR'
### Aliases: plot.bayesQR

### ** Examples

# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Analyze 5 quantiles using default prior
out <- bayesQR(y ~ X, quantile=c(.05,.25,.5,.75,.95), ndraw=5000)

# Check traceplot of first variable of .75 quantile regression 
plot(out, var=1, quantile=.75, plottype="trace")

# Check posterior histogram of first variable of .5 quantile regression 
plot(out, var=1, quantile=.5, plottype="hist")

# Create default quantile plot of first variable
plot(out, var=1, plottype="quantile")

# Create quantile plot of second variable with 90% credible interval
plot(out, var="X", credint=c(.05, .95), plottype="quantile", main="This is an example")



cleanEx()
nameEx("predict.bayesQR")
### * predict.bayesQR

flush(stderr()); flush(stdout())

### Name: predict.bayesQR
### Title: Calculate predicted probabilities for binary quantile regression
### Aliases: predict.bayesQR

### ** Examples

# Simulate data from binary regression model
set.seed(1234)
n <- 200
X <- matrix(runif(n=n*2, min=-5, max=5),ncol=2)
ystar <- cbind(1,X)%*% c(1,1.5,-.5) + rnorm(n=n, mean=0, sd=abs(2*X[,1]))
y <- as.numeric(ystar>0)

## Estimate a sequence of binary quantile regression models
out <- bayesQR(y ~ X, quantile=seq(.1,.9,.1), ndraw=4000)

# Calculate predicted probabilities
pred <- predict(object=out, X=cbind(1,X), burnin=2000)

# Make histogram of predicted probabilities 
hist(pred,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(pred>.5))

# Compare with logit model
mylogit <- glm(y ~ X, family=binomial(logit))

# Make histogram of predicted probabilities 
hist(mylogit$fit,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(mylogit$fit>.5))



cleanEx()
nameEx("print.bayesQR")
### * print.bayesQR

flush(stderr()); flush(stdout())

### Name: print.bayesQR.summary
### Title: Prints the contents of 'bayesQR.summary' object to the console
### Aliases: print.bayesQR.summary

### ** Examples

# Simulate data from heteroskedastic regression
set.seed(666)
n <- 200
X <- runif(n=n,min=0,max=10)
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Analyze 0.5 quantile using default prior and adaptive lasso
out <- bayesQR(y ~ X, alasso=TRUE, ndraw=5000) 

# Return Bayes estimates and credible intervals 
sum <- summary(out, burnin=1000)

# Print the bayesQR.summary object
sum



cleanEx()
nameEx("prior")
### * prior

flush(stderr()); flush(stdout())

### Name: prior
### Title: Create prior for Bayesian quantile regression
### Aliases: prior

### ** Examples

# Load the Prostate cancer dataset
data(Prostate)

# Create informative prior object
prior <- prior(lpsa~., data=Prostate, beta0=rep(5,9), V0=diag(9)) 

# Investigate structure of bayesQR.prior object
str(prior)

# Estimate the model parameters with informative prior
out <- bayesQR(lpsa~., data=Prostate, prior=prior, ndraw=5000)

# Print results
summary(out)



cleanEx()
nameEx("summary.bayesQR")
### * summary.bayesQR

flush(stderr()); flush(stdout())

### Name: summary.bayesQR
### Title: Summarize the output of the 'bayesQR' function
### Aliases: summary.bayesQR

### ** Examples

# Load the Prostate cancer dataset
data(Churn)

# Estimate the model parameters with default prior
out <- bayesQR(churn~gender+recency, data=Churn, ndraw=5000)

# Return Bayes estimates and credible intervals 
sum <- summary(out, burnin=1000)

# Inspect structure of bayesQR.summary object
str(sum)

# Print bayesQR.summary object
sum



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
