CalcProductLimitEst2 <- function(ageLast = ageLast, ageFirst = ageFirst,
                            departType = departType) {
  # Product limit estimator:
  idsort <- sort.int(ageLast, index.return = TRUE)$ix
  agev <- ageLast[idsort]
  nage <- length(agev)
  Cx <- rep(0, nage)
  delx <- rep(0, nage)
  for (ii in 1:nage) {
    idNx <- which(ageFirst <= agev[ii] & ageLast >= agev[ii])
    Cx[ii] <- length(idNx) / nage
    if (departType[idsort[ii]] == "D") delx[ii] <- 1
  }
  ple <- cumprod((1 - 1 / (nage * Cx))^delx)
  pleTab <- data.frame(Ages = agev, ple = ple)
  class(pleTab) <- "pointLimitEst"
  return(pleTab)
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


