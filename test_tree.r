###############################################################################
# A small experiment
library(huge)


###############################################################################
# TEST TCA
set.seed(10)
n = 2000
test_size = 500
d <- 3
a <- rep(c(3, 3, 3), d)

f <- function(x, alpha) {
  return(sign(x) * abs(x)^alpha)
}


L_true = huge.generator(n = n, d = d, graph="scale-free", v = 0.0)
s_true_full = L_true$data
#A_true = svd(matrix(rnorm(d^2), d, d))$u
A_true = diag(1, d)
W_true = solve(A_true)


# Regardless of which type of data was chosen, now transform it via the npn:
for(i in 1:d) {
  s_true_full[,i] <- f(s_true_full[,i], a[i]*0.2)
}



# Mix the data via the mixing matrix
ind <- sample(nrow(s_true_full), test_size)
s_true_train = s_true_full[-ind, ]
s_true_test = s_true_full[ind, ]
x_train = s_true_train %*% t(A_true)
x_test = s_true_test %*% t(A_true)

# Save x_train
write.table(t(x_train), "tca/x_train.csv", row.names=FALSE, col.names=FALSE, sep=",")

# Then run the TCA matlab code on it, save W, and read it in:
W_est = as.matrix(read.table("tca/W.csv", sep=","))


W_est %*% t(W_est)
