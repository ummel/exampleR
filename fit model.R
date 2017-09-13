# Load a package (this must be installed manually if packrat is enabled)
# This script need not be part of the published package
library(gbm)

# Create input dataset using the ?gbm example
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]
SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)
X1[sample(1:N,size=500)] <- NA
X4[sample(1:N,size=300)] <- NA
gbm.data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# Fit lm model to mtcars dataset
fitted_model <- gbm(Y ~ X1 + X2, data = gbm.data, keep.data = FALSE)

# Save model object to disk as .rda object
save(fitted_model, file = "data/fitted_model.rda")
