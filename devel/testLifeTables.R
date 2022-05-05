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

# theta <- c(-2, 3, 0.25, -8, 0.36)
# Calculate survival from selected model:
x <- seq(0, 100, 0.01)
Sx <- CalcSurv(theta = theta, x = x, model = model, shape = shape)

# number of individuals:
n <- 1000

# Study window:
timeStart <- 1970
studySpan <- 10 # 20
timeEnd <- timeStart + studySpan

# Simulate lifespans:
ageDeath <- SampleRandAge(n * 2, theta = theta, dx = 1 / 365.25,
                          model = model, shape = shape)

# ID of truncation:
# indTrunc <- rbinom(n = n, size = 1, 0.25)
# nTrunc <- sum(indTrunc)
# idTrunc <- which(indTrunc == 1)

# Age at truncation:
# ageFirst <- rep(0, n)
# ageFirst <- runif(n = nTrunc, min = rep(0, nTrunc),
#                   max = ageLast[idTrunc] / 2)

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

# Calculate life table:
lt <- CalcLifeTable(ageLast = ageLast, ageFirst = ageFirst, 
                    departType = departType)

# Calculate life table CIs:
# ltCIs <- CalcLifeTableCIs(ageLast = ageLast, ageFirst = ageFirst, 
#                     departType = departType)

# Plot life table CIs:
# plot(ltCIs, cex.lab = 1.25, cex.axis = 1.25, mar = c(2, 5, 1, 1),
#      cex.legend = 1.5, demorate = "all", lwd = 2)

# Calculate Kaplan-Meier curve:
KM <- CalcKaplanMeier(ageLast = ageLast, 
                      departType = departType)

# Kaplan-Meier CIs:
# KMcis <- CalcKaplanMeierCIs(ageLast = ageLast,
#                             departType = departType)

ple <- CalcProductLimitEst(ageLast = ageLast, ageFirst = ageFirst, 
                         departType = departType)
pleCIs <- CalcProductLimitEstCIs(ageLast = ageLast, ageFirst = ageFirst, 
                               departType = departType)

# Plot K-M and lx curves:
plot(lt[, c("Ages", "lx")], type = 'p', col = 'red', lwd = 4, ylim = c(0, 1),
     pch = 19)
# polygon(c(ltCIs$lx[, "Ages"], rev(ltCIs$lx[, "Ages"])), 
#         c(ltCIs$lx[, 'Lower'], rev(ltCIs$lx[, "Upper"])), 
#         col = adjustcolor(col = 'red', alpha.f = 0.25), border = NA)
# 
lines(KM$Ages, KM$KM, col = 'orange', type = 's', lwd = 2)
# polygon(c(KMcis$KM$Ages, rev(KMcis$KM$Ages)), 
#         c(KMcis$KM$Lower, rev(KMcis$KM$Upper)), 
#         col = adjustcolor(col = 'orange', alpha.f = 0.25), border = NA)

lines(x, Sx, col = 'dark red', lwd = 2)

# plot(lt, demorate = "lx", cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
# 
# plot(lt, demorate = "qx", cex.axis = 1.25, cex.lab = 1.25, lwd = 2)

CalcAgeingRate(x = c(10, 20), theta = theta, model = model, shape = shape)

ex <- CalcRemainLifeExp(x = c(0:ceiling(max(ageLast))), theta = theta, model = model, shape = shape, 
                  dx = 0.001)

xvv <- seq(0, ceiling(max(ageLast)), 0.01)
# Plot mortality:

par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))
# Plot mortality:
plot(xvv, CalcMort(theta = theta, x = xvv, model = model, shape = shape),
     type = 'l', xlab = "", ylab = "Hazard rate", col = 2, lwd = 2)

# Survival:
plot(c(0, max(ageLast)), c(0, 1), col = NA, xlab = "", ylab = "Survival")
lines(xvv, CalcSurv(theta = theta, x = xvv, model = model, shape = shape), col = 'dark green', lwd = 2)
xlt <- ltCIs$lx[, "Ages"]
ylt <- ltCIs$lx[, c("Lower", "Upper")]
polygon(c(xlt, rev(xlt)), c(ylt[, 1], rev(ylt[, 2])), 
        col = adjustcolor('red', alpha.f = 0.25), border = NA)
lines(lt[, c("Ages", 'lx')], type = 's', col = 'red')
lines(kmTab$Age, kmTab$KM, col = 3, type = 's')

# Remaining life expectancy:
plot(lt[, c("Ages", "ex")], type = 'l')
xci <- ltCIs$ex[, "Ages"]
yci <- ltCIs$ex[, c("Lower", "Upper")]
polygon(c(xci, rev(xci)), c(yci[, 1], rev(yci[, 2])), col = adjustcolor('grey40', alpha.f = 0.25), border = NA)
lines(ex, col = 2)
