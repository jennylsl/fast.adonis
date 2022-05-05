library(vegan)
library(plyr)
library(fast.adonis)
set.seed(10010)


# no weights
# create a distance matrix
dim <- 4
D <- dist(matrix(rnorm(10^dim*2), ncol=10))
# transform to A
A <- -0.5*as.matrix(D)^2
# dim(A)
# D <- NULL
# create a matrix for independent variables
X <- data.frame(matrix(rnorm(10^dim*2), ncol=10))
# fast adonis
fit0 <- fast.adonis(A ~ X1, data=X, permutations = 50, boot.times = 10)
fit1 <- fast.adonis(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=X, permutations = 50, boot.times = 10)
# fit0$aov.tab
# fit1$aov.tab

# with weights
weights <- abs(rnorm(dim(A)[1]))
fit2 <- fast.adonis(A ~ X1, data=X, permutations = 50, boot.times = 10, weights = weights)
# fit2$aov.tab


