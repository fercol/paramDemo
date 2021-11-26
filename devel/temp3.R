
# Calculate age at maximum fecundity:
CalcAgeMaxFec <- function(beta, modelFec = "M1") {

  # a) Fecundity method:
  .CalcFecund <- function(beta, ...) UseMethod(".CalcFecund")
  .CalcFecund.matrix <- .DefineFecundityMatrix(modelFec = modelFec)
  .CalcFecund.numeric <- .DefineFecundityNumeric(modelFec = modelFec)
  
  # Find age at maximum fecundity:
  xv <- seq(0, 100, 0.0001)
  if (modelFec == "M1") {
    xm <- beta["b2"]
    dd <- 0
    ii <- 0
  } else if (modelFec %in% c("M2", "M3")) {
    if (modelFec == "M2") {
      dfdx <- function(x, beta) {
        w <- x - beta["b2"]
        u <- exp(-beta["b1"] * w^2)
        z <- exp(-beta["b3"] * w)
        v <- 1 + z
        du <- -2 * beta["b1"] * w * u
        dv <- -beta["b3"] * z
        df <- (du * v - u * dv) / v^2
        return(df)
      }
      dfdx2 <- function(x, beta) {
        w <- x - beta["b2"]
        u <- exp(-beta["b1"] * w^2)
        z <- exp(-beta["b3"] * w)
        v <- 1 + z
        du <- -2 * beta["b1"] * w * u
        dv <- -beta["b3"] * z
        ddu <- (-2 * beta["b1"] + (-2 * beta["b1"] * w)^2) * u
        ddv <- beta["b3"]^2 * z
        df <- (du * v - u * dv) / v^2
        df2 <- (ddu * v - u * ddv) / v^2 - 2 * dv / v * df
        return(df2)
      }
    } else {
      dfdx <- function(x0, beta) {
        x0^3 + x0^2 * (2 - beta["b2"]) + x0 * (1 - 2 * beta["b2"]) +
          beta["b3"] / (2 * beta["b1"]) - beta["b2"]
      }
      dfdx2 <- function(x0, beta) {
        3 * x0^2 + 2 * x0 * (2 - beta["b2"]) + 1 - 2 * beta["b2"]
      }
    }
    id0 <- which(sign(dfdx(xv[-length(xv)], beta)) != sign(dfdx(xv[-1], beta)))
    if (length(id0) > 1) {
      id0 <- id0[which(abs(xv[id0] - beta["b2"]) == 
                         min(abs(xv[id0] - beta["b2"])))]
    }
    x0 <- xv[id0]
    dd <- 1
    ii <- 0
    while (dd > 0.0001 & ii <= 25) {
      ii <- ii + 1
      x1 <- x0 - dfdx(x0, beta) / dfdx2(x0, beta)
      dd <- abs(x1 - x0)
      x0 <- x1
    }
    if (dd < 0.0001) {
      xm <- x1
    } else {
      xm <- NA
    }
  }
  maxFec <- .CalcFecund(beta, xm)
  maxFecv <- c(xm, maxFec, dd, ii)
  names(maxFecv) <- c("Age", "Fecund", "error", "itterations")
  return(maxFecv)
}


xxv <- seq(0, 20, 0.001)
beta <- c(b0 = 0.75, b1 = 0.025, b2 = 10, b3 = -20)
fec <- CalcFecund(beta, xxv, modelFec = "M3", checkBeta = T)
maxfec <- CalcAgeMaxFec(beta = beta, modelFec = "M3")
plot(xxv, fec, type = 'l')
abline(v = maxfec["Age"], col = 2, lty = 2)
abline(h = maxfec["Fecund"], col = 2, lty = 2)


beta <- c(b0 = 0.75, b1 = 0.005, b2 = 8, b3 = -0.2)
fec <- CalcFecund(beta, xxv, modelFec = "M2", checkBeta = T)
maxfec <- CalcAgeMaxFec(beta = beta, modelFec = "M2")
plot(xxv, fec, type = 'l')
abline(v = maxfec["Age"], col = 2, lty = 2)
abline(h = maxfec["Fecund"], col = 2, lty = 2)

fecund <- function(beta, x) {
  fec <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2 + 
                            beta["b3"] * 1/(x + beta["b4"]))
  return(fec)
}


beta <- c(b0 = 0.775, b1 = 0.002, b2 = 4, b3 = -0.5)
fec <- CalcFecund(beta, xxv, modelFec = "M3", checkBeta = T)

#betatemp <- c(beta, b4 = 5)
fectemp <- fecund(betatemp, xxv)
plot(xxv, fec, type = 'l')
lines(xxv, fectemp, col = 2)


# MORTALITY RATE DOUBLING TIME (MRDT):
.CalcMRDT <- function(theta, x, model = "GO", shape = "simple",
                      checkTheta = TRUE, error = 1e-10) {
  
  # Verify model and shape:
  .VerifyModelShape(model = model, shape = shape)
  
  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }
  
  # Calculate MRDT:
  if (shape == "simple") {
    if (model == "EX") {
      MRDT <- 0
      ageDep <- FALSE
    } else if (model == "GO") {
      MRDT <- log(2) / theta["b1"]
      ageDep <- FALSE
    } else if (model == "WE") {
      MRDT <- x * (2^(1/(theta["b0"] - 1)) - 1)
      ageDep <- TRUE
    } else {
      g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
      M <- (2 - 2 * g) / (1 - g * (exp(theta["b1"] * x) + 1))
      MRDT <- log(M) / theta["b1"]
      ageDep <- TRUE
    }
    errMRDT <- 0
    iter <- 0
  } else if (shape == "Makeham") {
    if (model == "GO") {
      MRDT <- log(theta["c"] / (exp(theta["b0"] + theta["b1"] * x)) + 2) /
        theta["b1"]
      ageDep <- TRUE
    } else if (model == "WE") {
      MRDT <- (theta["c"] / (theta["b0"] * theta["b1"]^(theta["b0"])) + 
                 2 * x^(theta["b0"] - 1))^(1/(theta["b0"] - 1)) - x
      ageDep <- TRUE
    } else {
      g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
      M <- theta["c"] / exp(theta["b0"] + theta["b1"] * x) + 
        2 / (1 + g * (exp(theta["b1"] * x) - 1))
      MM <- (M * (1 - g)) / (1 - M * g * exp(theta["b1"] * x))
      MRDT <- log(abs(MM)) / theta["b1"]
      ageDep <- TRUE
    }
    errMRDT <- 0
    iter <- 0
  } else {
    xv <- seq(0, 1000, 0.01)
    mux <- CalcMort(theta = theta, x = x, model = model, shape = shape)
    errMRDT <- 1
    iter <- 0
    while(errMRDT > error & iter <= 25) {
      iter <- iter + 1
      muv <- CalcMort(theta = theta, x = xv, model = model, shape = shape)
      id <- which(abs(muv - 2 * mux) == min(abs(muv - 2 * mux)))
      errMRDT <- abs(muv[id] - 2 * mux)
      xmrdt <- xv[id]
      xv <- seq(xv[id-1], xv[id+1], length = 1000)
    }
    MRDT <- xmrdt - x
  }
  
  MRDTout <- c(x, MRDT, errMRDT, iter)
  names(MRDTout) <- c("IniAge", "MRDT", "Error", "Iterations")
  return(MRDTout)
}


TestMRDTfun <- function(theta, x, model = "GO", shape = "simple") {
  MRDT <- .CalcMRDT(theta = theta, x = x, model = model, shape = shape)
  mux <- CalcMort(theta = theta, x = x, model = model, shape = shape)
  muxh <- CalcMort(theta = theta, x = x + MRDT["MRDT"], model = model, 
                   shape = shape)
  cat(sprintf("mu ratio = %s\n", round(muxh / mux, 3)))
  cat(sprintf("MRDT     = %s years\n", round(MRDT["MRDT"], 3)))
}

x <- 5
# GO-simple:
theta <- c(b0 = -5, b1 = 0.15)
dem <- CalcDemo(theta = theta)
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta, x)

# LO-simple:
x <- 50
theta <- c(b0 = -10, b1 = 0.15, b2 = 2)
dem <- CalcDemo(theta = theta, model = "LO")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "LO")
xv <- seq(0, max(dem$age), length = 10000)
g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
M <- (2 - 2 * g) / (1 - g * (exp(theta["b1"] * xv) + 1))
vazym <- xv[which(sign(M[-1])!=sign(M[-length()]))]
plot(xv, M)

# GO-Makeham:
theta <- c(c = 0.001, b0 = -3, b1 = 0.15)
dem <- CalcDemo(theta = theta, shape = "Makeham")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, shape = "Makeham")

# WE-Makeham:
theta <- c(c = 0.1, b0 = 1.5, b1 = 0.2)
dem <- CalcDemo(theta = theta, model = "WE", shape = "Makeham")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "WE", shape = "Makeham")

# LO-Makeham:
x <- 50
theta <- c(c = 0.001, b0 = -10, b1 = 0.15, b2 = 2)
dem <- CalcDemo(theta = theta, model = "LO", shape = "Makeham")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "LO", shape = "Makeham")

# GO-bathtub:
theta <- c(a0 = 0, a1 = 2, c = 0.001, b0 = -5, b1 = 0.15)
dem <- CalcDemo(theta = theta, model = "GO", shape = "bathtub")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "GO", shape = "bathtub")

# WE-bathtub:
theta <- c(a0 = 0, a1 = 2, c = 0.001, b0 = 1.5, b1 = 0.2)
dem <- CalcDemo(theta = theta, model = "WE", shape = "bathtub")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "WE", shape = "bathtub")

# LO-bathtub:
theta <- c(a0 = 1, a1 = 2, c = 0.001, b0 = -5, b1 = 0.15, b2 = 2)
dem <- CalcDemo(theta = theta, model = "LO", shape = "bathtub")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "LO", shape = "bathtub")

# ==== test code: ====
xv <- seq(0, 1000, 0.01)
theta <- c(a0 = 1, a1 = 2, c = 0.001, b0 = -5, b1 = 0.15, b2 = 2)
mux <- CalcMort(theta = theta, x = x, model = "LO", shape = 'bathtub')
err <- 1
iter <- 0
while(err > 0.0000001) {
  iter <- iter + 1
  muv <- CalcMort(theta = theta, x = xv, model = "LO", shape = 'bathtub')
  id <- which(abs(muv - 2 * mux) == min(abs(muv - 2 * mux)))
  err <- abs(muv[id] - 2 * mux)
  xmrdt <- xv[id]
  xv <- seq(xv[id-1], xv[id+1], length = 1000)
}

muxh <- CalcMort(theta = theta, x = xmrdt, model = "LO", 
                 shape = 'bathtub')
muxh / mux

theta <- c(c = 0.01, b0 = -5, b1 = 0.15, b2 = 2)

g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
M <- theta["c"] / exp(theta["b0"] + theta["b1"] * x) + 
  2 / (1 + g * (exp(theta["b1"] * x) - 1))
MM <- (M * (1 - g)) / (1 - M * g * exp(theta["b1"] * x))
MRDT <- log(MM) / theta["b1"]



dmudx <- function(theta, x) {
  al <- exp(theta["b0"])
  g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
  b1 <- theta["b1"]
  eb1x <- exp(b1 * x)
  denom <- 1 + g * (eb1x - 1)
  dmu <- al * b1 * eb1x * (1 - g) / denom^2
  return(dmu)
}

dmudx2 <- function(theta, x) {
  al <- exp(theta["b0"])
  g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
  b1 <- theta["b1"]
  eb1x <- exp(b1 * x)
  denom <- 1 + g * (eb1x - 1)
  dmu2 <- al * b1^2 * eb1x * (1-g) * denom * (1 - g * (eb1x + 1)) / denom^4
  return(dmu2)
}

theta <- c(b0 = -3, b1 = 0.15, b2 = 0.2)
xv <- seq(0, 50, 0.01)
g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
b1 <- theta["b1"]
xinfl2 <- log(abs(g^2 + 2 * g - 1) / g^2) / (2 * b1)
dmu2 <- dmudx2(theta, xv)
xinfl <- xv[which(sign(dmu2[-1]) != sign(dmu2[-length(dmu2)]))]

par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))
plot(xv, CalcMort(theta, x = xv, model = "LO"), type = 'l', lwd = 4, 
     col = 'dark red')
abline(v = xinfl, col = 'orange', lty = 3, lwd = 6)
abline(v = xinfl2, col = 'red', lty = 2, lwd = 2)
plot(xv, dmudx(theta, xv), type = 'l', lwd = 4, col = 'dark red')
abline(v = xinfl, col = 'orange', lty = 3, lwd = 6)
abline(v = xinfl2, col = 'red', lty = 2, lwd = 2)
plot(xv, dmudx2(theta, xv), type = 'l', lwd = 4, col = 'dark red')
abline(v = xinfl, col = 'orange', lty = 3, lwd = 6)
abline(v = xinfl2, col = 'red', lty = 2, lwd = 2)


# LO-simple:
x <- 17.623759
theta <- c(b0 = -3, b1 = 0.05, b2 = 2)
dem <- CalcDemo(theta = theta, model = "LO")
plot(dem$age, dem$mort, type = 'l', col = 'dark red', lwd = 4)
abline(v = x, lty = 2, col = 'orange', lwd = 4)
TestMRDTfun(theta = theta, x = x, model = "LO")
xv <- seq(0, max(dem$age), length = 10000)
b1 <- theta["b1"]
xinfl2 <- log(abs(g^2 + 2 * g - 1) / g^2) / (2 * b1)

g <- theta["b2"] * exp(theta["b0"]) / theta["b1"]
M <- (2 - 2 * g) / (1 - g * (exp(b1 * xv) + 1))
if (g > 1 / (1 + exp(b1 * x)) & g < 1) {
  xazym <- log((1-g)/g) / b1
} else {g > 1}


vazym <- xv[which(sign(M[-1])!=sign(M[-length(M)]))]
plot(xv, M)
abline(v = vazym, col = "orange")
abline(v = xinfl2, col = 'red')

xazym <- log((1-g)/g) / b1
g > 1 / (1 + exp(b1 * x))
g < 1
