CalcCIs <- function(n, ns, alpha = 0.05, method = "ClopperPearson") {
  ln <- length(ns)
  if (method == "ClopperPearson" | method == "clopperpearson") {
    Lower <- qbeta(alpha / 2, ns, n - ns + 1)
    Upper <- qbeta(1 - alpha / 2, ns + 1, n - ns)
    cis <- cbind(Lower, Upper)
  } else {
    z <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
    cis <- (ns + z^2 / 2) / (n + z^2) + z / (n + z^2) * 
      sqrt((ns * (n - ns)) / n + z^2 / 4) * 
      matrix(c(-1, 1), ln, 2, byrow = TRUE)
    colnames(cis) <- c("Lower", "Upper")
    cis[which(cis[, "Lower"] < 0), "Lower"] <- 0
    cis[which(cis[, "Upper"] > 1), "Upper"] <- 1
  }
  return(cis)
}
