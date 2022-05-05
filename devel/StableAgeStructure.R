# Function to construct Leslie matrix:
FillMatr <- function(p, f, n) {
  idcol <- 1:n
  idrow <- c(2:n, n)
  idpa <- (idcol - 1) * n + idrow
  Aa <- matrix(0, n, n)
  Aa[1, ] <- f
  Aa[idpa] <- p
  return(Aa)
}
# Using Euler-Lotka's renewal equation:
Findr <- function(r) {
  lhs <- sum(exp(-r * xv) * Sxv * mxv * dx)
  return((lhs - 1)^2)
}

# Parameters:
theta <- c(-1, 1, 0.001, -6, 0.15)
beta <- c(0.2, 0.002, 20)

dx <- 0.0001
xv <- seq(0, 50, dx)
Sxx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
maxX <- ceiling(x[which(Sxx <= 0.0001)[1]])
x <- seq(0, maxX, dx)
nx <- length(x)
alpha <- 5
Sxv <- CalcSurv(theta = theta, x = xv, model = "GO", shape = "bathtub")
mxv <- c(rep(0, length(which(xv < alpha))), 
            CalcFert(beta = beta, x = xv[xv >= alpha], modelFert = "M1"))

rout <- optimize(f = Findr, interval = c(0, 5), tol = 10e-20)
r1 <- rout$minimum

x <- 0:maxX
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
mx <- c(rep(0, length(which(x < alpha))), 
        CalcFert(beta = beta, x = x[x >= alpha], modelFert = "M1"))
plot(x, Sx)
plot(x, mx)
px <- CalcSurv(theta = theta, x = x + 1, model = "GO", shape = "bathtub") /
  CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
w1 <- (Sx * exp(-r1 * x)) / sum(Sx * exp(-r1 * x))

# Using matrix algebra:
A <- FillMatr(p = px, f = mx, n = length(px))
eA <- eigen(A)
lambda <- Re(eA$value[1])
w2 <- abs(Re(eA$vector[, 1])) / sum(abs(Re(eA$vector[, 1])))
r2 <- log(lambda)

# Matrix with larger number of ages (sub-ages):
dx2 <- 0.5
xv2 <- seq(0, maxX, 0.5)
Sx2 <- CalcSurv(theta = theta, x = xv2, model = "GO", shape = "bathtub")
mx2 <- c(rep(0, length(which(xv2 < alpha))), 
        CalcFert(beta = beta, x = xv2[xv2 >= alpha], modelFert = "M1"))
px2 <- CalcSurv(theta = theta, x = xv2 + dx2, model = "GO", shape = "bathtub") /
  CalcSurv(theta = theta, x = xv2, model = "GO", shape = "bathtub")

A2 <- FillMatr(px2, mx2, length(px2))
eA2 <- eigen(A2)
lambda2 <- Re(eA2$value[1])
r3 <- log(lambda2)

# Combining both:
w3 <- (Sx * exp(-r2 * x)) / sum(Sx * exp(-r2 * x))

# Compare the age structures:
par(mfrow = c(2, 1))
plot(x, px, pch = 19, ylim = c(0, 1))
points(x[x>= alpha], mx[x >= alpha], col = 2, pch = 15)
plot(range(x), c(0, max(c(w1, w2))), col = NA, xlab = "Age", 
     ylab = "Stable age struct.")
lines(x, w1, col = 2, lwd = 4)
lines(x, w2, col = 3, lwd = 2)
lines(x, w3, col = 4, lwd = 1)
r1
r2
# Project the population and find r:
N <- rep(1, nx)
nt <- 1000
M <- rep(0, nt)
for (tt in 1:nt) {
  N <- A %*% N
  M[tt] <- sum(N)
}
lamb2 <- c(M[-1] / M[-nt])[nt-1]
log(lamb2)
