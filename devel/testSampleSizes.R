source("~/FERNANDO/PROJECTS/4.PACKAGES/paramDemo/pkg/R/paramDemo.R")

scen <- "GB"
if (scen == "GS") {
  # Gompertz:
  theta <- c(b0 = -5, b1 = 0.1)
  model <- "GO"
  shape <- "simple"
} else if (scen == "GM") {
  # Gompertz-Makeham:
  theta <- c(c = 0.01, b0 = -5, b1 = 0.1)
  model <- "GO"
  shape <- "Makeham"
} else if (scen == "GB") {
  # Siler:
  theta <- c(a0 = 0, a1 = 1, c = 0.001, b0 = -5, 
             b1 = 0.1)
  model <- "GO"
  shape <- "bathtub"
} else if (scen == "WS") {
  theta <- c(b0 = 2.5, b1 = 0.1)
  model <- "WE"
  shape <- "simple"
} else if (scen == "WM") {
  theta <- c(c = 0.001, b0 = 2.5, b1 = 0.1)
  model <- "WE"
  shape <- "Makeham"
} else if (scen == "WB") {
  theta <- c(a0 = -1, a1 = 2, c = 0.001, b0 = 2.5, b1 = 0.1)
  model <- "WE"
  shape <- "bathtub"
} else if (scen == "LS") {
  # Gompertz:
  theta <- c(b0 = -5, b1 = 0.1, b2 = 0.5)
  model <- "LO"
  shape <- "simple"
} else if (scen == "LM") {
  # Gompertz-Makeham:
  theta <- c(c = 0.01, b0 = -5, b1 = 0.1, b2 = 0.5)
  model <- "LO"
  shape <- "Makeham"
} else if (scen == "LB") {
  # Siler:
  theta <- c(a0 = -2, a1 = 1, c = 0.001, b0 = -10, 
             b1 = 0.1, b2 = 0.5)
  model <- "LO"
  shape <- "bathtub"
}

# Calculate survival from selected model:
x <- seq(0, 100, 0.01)
Sx <- CalcSurv(theta = theta, x = x, model = model, shape = shape)

# number of individuals:
n <- 30

# Study window:
timeStart <- 1970
studySpan <- 20
timeEnd <- timeStart + studySpan

# Number of iterations:
niter <- 1000

# Output matrix:
outmat <- matrix(0, niter, 3, dimnames = list(NULL, c("RMSE", "PoutCIs", "N")))

# Progress bar:
pb <- txtProgressBar(min = 1, max = niter, initial = 1)

# Start loop:
Start <- Sys.time()
for (iter in 1:niter) {
  # Simulate lifespans:
  ageDeath <- SampleRandAge(n * 3, theta = theta, dx = 1 / 365.25,
                            model = model, shape = shape)
  
  # Simulate times of birth:
  birthDate <- runif(n = n * 3, min = -studySpan, 
                     max = studySpan)
  
  # Times of death:
  deathDate <- birthDate + ageDeath
  
  # Find individuals that were alive after study started:
  idincl <- which(deathDate >= 0 & birthDate < studySpan)
  N <- length(idincl)
  
  # Reduce N to sample sizes:
  if (N != n) {
    idincl <- sample(idincl, n, replace = FALSE)
  }
  N <- length(idincl)
  
  # Subset dataset to detected individuals:
  ageDeath <- ageDeath[idincl]
  birthDate <- birthDate[idincl]
  deathDate <- deathDate[idincl]
  
  # Assign age first:
  ageFirst <- rep(0, N)
  ageFirst[birthDate < 0] <- abs(birthDate[birthDate < 0])
  
  # Find censored individuals:
  idcens <- which(deathDate > studySpan)
  
  # Assign age last:
  ageLast <- ageDeath
  ageLast[idcens] <- (studySpan - birthDate[idcens])
  
  # Assign depart type:
  departType <- rep("D", N)
  departType[idcens] <- "C"
  
  # Calculate Kaplan-Meier curve:
  KM <- CalcKaplanMeier(ageLast = ageLast, ageFirst = ageFirst, 
                        departType = departType)
  
  # Kaplan-Meier CIs:
  KMcis <- CalcKaplanMeierCIs(ageLast = ageLast, ageFirst = ageFirst, 
                              departType = departType)
  
  # Number of ages:
  nAges <- nrow(KMcis$KM)
  
  # Calculate Survival at KM ages:
  Sx <- CalcSurv(theta = theta, x = KMcis$KM$Ages, model = model, shape = shape)
  
  # Calculate root mean square errors:
  RMSE <- sqrt(sum((Sx - KMcis$KM$KM)^2) / N)
  
  # Find proportion of Sx that fall out of the K-M CIs:
  PoutCIs <- length(which(Sx < KMcis$KM$Lower | Sx > KMcis$KM$Upper)) / nAges
  
  # store results:
  outmat[iter, ] <- c(RMSE, PoutCIs, N)
  
  # Progress bar:
  setTxtProgressBar(pb = pb, value = iter)
}
End <- Sys.time()
End - Start

# Extract results:
QuanRes <- apply(outmat[, c("RMSE", "PoutCIs")], 2, quantile, 
                 c(0.025, 0.975), na.rm = TRUE)
# Plot results:
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(density(outmat[, "RMSE"], na.rm = TRUE), main = "RMSE", 
     xlim = QuanRes[, "RMSE"])
plot(density(outmat[, "PoutCIs"]), main = "PoutCIs",
     xlim = QuanRes[, "PoutCIs"])

