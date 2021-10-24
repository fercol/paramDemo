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

# number of individuals:
n <- 10000

# Simulate lifespans:
ageLast <- SampleRandAge(n, theta = theta, dx = 1 / 365.25,
                         model = model, shape = shape)
ageFirst <- rep(0, n)
departType <- rep("D", n)

lt <- CalcLifeTable(ageLast = ageLast, ageFirst = ageFirst, 
                    departType = departType)

ltCIs <- CalcLifeTableCIs(ageLast = ageLast, ageFirst = ageFirst, 
                    departType = departType)


plot(ltCIs, cex.lab = 1.25, cex.axis = 1.25, mar = c(2, 5, 1, 1),
     cex.legend = 1.5, demorate = "all", lwd = 2)


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

# Remaining life expectancy:
plot(lt[, c("Ages", "ex")], type = 'l')
xci <- ltCIs$ex[, "Ages"]
yci <- ltCIs$ex[, c("Lower", "Upper")]
polygon(c(xci, rev(xci)), c(yci[, 1], rev(yci[, 2])), col = adjustcolor('grey40', alpha.f = 0.25), border = NA)
lines(ex, col = 2)
