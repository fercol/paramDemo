library(survival)
source("~/FERNANDO/PROJECTS/4.PACKAGES/paramDemo/pkg/R/paramDemo.R")
n <- 1000
x <- SampleRandAge(n, theta = c(b0 = -3, b1 = 0.1), dx = 0.0001, model = 'GO',
                   shape = 'simple')
ev <- rep(1, n)
tab <- data.frame(ID = 1:n, Age = x, ev = ev)
Sxcalc <- Surv(time = rep(0, n), time2 = x, event = ev)
fit <- survfit(Sxcalc ~ 1, data = tab, id = tab$ID)
plot(fit)
Sx <- exp(-fit$cumhaz)
ages <- fit$time
ex <- sum(diff(ages) * Sx[-length(Sx)])
lt <- CalcLifeTable(ageLast = x, ageFirst = rep(0, n), departType = rep("D", n))


alpha <- 2
x2 <- x[which(x > alpha)]
n2 <- length(x2)
tab2 <- data.frame(ID = 1:n2, Age = x2, ev = rep(1, n2))
Sxcalc2 <- Surv(time = rep(0, n2), time2 = x2, event = tab2$ev)
fit <- survfit(Sxcalc2 ~ 1, data = tab2, id = tab2$ID)
plot(fit)

