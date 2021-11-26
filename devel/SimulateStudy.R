# ============================ CODE METADATA ================================= #
# PACKAGE: paramDemo (ancillary functions)
# AUTHOR: Fernando Colchero
# DATE CREATED: 2020-06-01
# DATE MODIFIED: 2021-10-14
# DESCRIPTION: Ancillary script to simulate a population in time.
# COMMENTS: 
# ============================ START CODE ==================================== #
source("~/FERNANDO/PROJECTS/4.PACKAGES/paramDemo/pkg/R/paramDemo.R")

scen <- "WB"
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
  theta <- c(a0 = -2, a1 = 1, c = 0.001, b0 = -10, 
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
n <- 1000

# Study window:
timeStart <- 1970
studySpan <- 20
timeEnd <- timeStart + studySpan

# Simulate lifespans:
ageDeath <- SampleRandAge(n * 2, theta = theta, dx = 1 / 365.25,
                          model = model, shape = shape)

# Simulate times of birth:
birthDate <- runif(n = n * 2, min = -studySpan, 
                   max = studySpan)

# Times of death:
deathDate <- birthDate + ageDeath

# Find individuals that were alive after study started:
idincl <- which(deathDate >= 0 & birthDate < studySpan)
n <- length(idincl)

# Subset dataset to detected individuals:
ageDeath <- ageDeath[idincl]
birthDate <- birthDate[idincl]
deathDate <- deathDate[idincl]

# Assign age first:
ageFirst <- rep(0, n)
ageFirst[birthDate < 0] <- abs(birthDate[birthDate < 0])

# Find censored individuals:
idcens <- which(deathDate > studySpan)

# Assign age last:
ageLast <- ageDeath
ageLast[idcens] <- (studySpan - birthDate[idcens])

# Assign depart type:
departType <- rep("D", n)
departType[idcens] <- "C"
