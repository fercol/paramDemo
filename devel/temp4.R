CalcKaplanMeier2 <- function(ageLast = ageLast, ageFirst = ageFirst,
                            departType = departType) {
  unAgeSort <- sort(unique(ageLast[which(departType == "D")]))
  nAges <- length(unAgeSort)
  kmRatio <- rep(0, nAges)
  for (aa in 1:nAges) {
    age <- unAgeSort[aa]
    nDead <- length(which(ageLast == age & departType == "D"))
    nTot <- length(which(ageLast >= age & ageFirst <= age))
    kmRatio[aa] <- nDead / nTot
  }
  km <- cumprod(1 - kmRatio)
  kmTab <- data.frame(Age = unAgeSort, KM = km)
  class(kmTab) <- "KaplanMeier"
  return(kmTab)
}

Start1 <- Sys.time()
out1 <- CalcKaplanMeier(ageLast = ageLast, ageFirst = ageFirst, 
                        departType = departType)
End1 <- Sys.time()

Start2 <- Sys.time()
out2 <- CalcKaplanMeier2(ageLast = ageLast, ageFirst = ageFirst, 
                        departType = departType)
End2 <- Sys.time()

End1 - Start1
End2 - Start2

par(mfrow = c(2, 1))
plot(out1$Age, out1$KM, col = 'orange', lwd = 4, type = 'l')
lines(out2$Age, out2$KM, col = 'dark red', lwd = 1)
plot(out1$KM - out2$KM, type= 'l')

# Calculation of life table CIs:
CalcLifeTableCIs <- function(ageLast, ageFirst = NULL, departType, 
                             nboot = 2000, alpha = 0.05) {
  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, n)
  }
  
  # Unit age vector for that sex:
  agev <- 0:ceiling(max(ageLast))
  nage <- length(agev)
  
  # Mean life table:
  ltMean <- CalcLifeTable(ageLast = ageLast, ageFirst = ageFirst, 
                          departType = departType)
  
  # Labels for demographic rates:
  demRates <- c("lx", "qx", "px", "ex")
  
  # Output bootstrap array
  bootarr <- array(0, dim = c(nage, nboot, 4), 
                   dimnames = list(NULL, NULL, demRates))
  
  # Run bootstrap:
  for (iboot in 1:nboot) {
    idboot <- sample(1:n, size = n, replace = TRUE)
    ageLastBoot <- ageLast[idboot]
    ageFirstBoot <- ageFirst[idboot]
    departTypeBoot <- departType[idboot]
    ltb <- CalcLifeTable(ageLast = ageLastBoot, ageFirst = ageFirstBoot, 
                         departType = departTypeBoot)
    nl <- nrow(ltb)
    for (dr in demRates) {
      bootarr[1:nl, iboot, dr] <- ltb[, dr]
    }
  }
  
  # Calculate CIs:
  ltcis <- sapply(demRates, function(dr) {
    cii <- t(apply(bootarr[, , dr], 1, quantile, c(alpha / 2, 1 - alpha / 2),
                   na.rm = TRUE))
    colnames(cii) <- c("Lower", "Upper")
    return(cii)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # Merge with mean life table:
  for (dr in demRates) {
    ltcis[[dr]] <- cbind(ltMean[, c("Ages", dr)], ltcis[[dr]])
  }
  
  # Settings:
  ltcis$Settings <- c(nboot = nboot, alpha = alpha)
  
  # class:
  class(ltcis) <- "lifeTableCIs"
  return(ltcis)
}

