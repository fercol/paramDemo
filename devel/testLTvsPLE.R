load("devel/testLTvsPLE.RData")
ageLast <- testdf$ageLast
ageFirst <- testdf$ageFirst
departType <- testdf$departType

# =========== #
# LIFE TABLE:
# =========== #
dx <- 0.5
# Number of records:
n <- length(ageLast)

# Set age first to 0 if NULL:
if (is.null(ageFirst)) {
  ageFirst <- rep(0, n)
}

# Unit age vector for that sex:
agev <- seq(from = 0, to = ceiling(max(ageLast[which(departType == "D")])), 
            by = dx)
nage <- length(agev)

# Outputs:
Nx <- Dx <- ax <- rep(0, nage)
for (xx in 1:nage) {
  # A) EXPOSURES:
  # Find how many entered the interval (including truncated):
  idNx <- which(ageFirst < agev[xx] + dx & ageLast >= agev[xx])
  
  # Extract ages and departType:
  xf <- ageFirst[idNx]
  xl <- ageLast[idNx]
  dt <- departType[idNx]
  
  # proportion of truncation in interval:
  trp <- xf - agev[xx]
  trp[trp < 0] <- 0
  
  # proportion of censoring:
  cep <- agev[xx] + dx - xl
  cep[cep < 0] <- 0
  cep[dt == "D"] <- 0
  
  # Calculate exposures:
  nexp <- 1 - trp - cep
  Nx[xx] <- sum(nexp)
  
  # B) DEATHS:
  # Calculate total deaths in the interval:
  idDx <- which(dt == "D" & xl < agev[xx] + dx)
  # Dx[xx] <- length(idDx)
  Dx[xx] <- sum(nexp[idDx])
  
  # C) PROPORTION LIVED BY THOSE THAT DIED IN INTERVAL:
  if (Dx[xx] > 1) {
    ylived <- xl[idDx] - agev[xx]
    ax[xx] <- sum(ylived) / Dx[xx]
  } else {
    ax[xx] <- 0
  }
}

# Age-specific mortality probability:
qx <- Dx / Nx
qx[which(is.na(qx))] <- 0

# Age-specific survival probability:
px <- 1 - qx

# Survivorship (or cumulative survival):
lx <- c(cumprod(px))[1:nage]
# lx <- Nx / n

# Number of individual years lived within the interval:
Lx <- lx * (1 - ax * qx)
# Note: correction on the calculation of Lx (doesn't work when
#       there are censored and truncated records)
# Lx <- Nx - Dx * ax
Lx[is.na(Lx)] <- 0

# Total number of individual years lived after age x:
Tx <- rev(cumsum(rev(Lx))) * dx

# Remaining life expectancy after age x:
ex <- Tx / lx 
ex[which(is.na(ex))] <- 0
# (Note: follows on the correction for Lx)
# ex <- Tx / Nx

# Life-table:
lt <- data.frame(Ages = agev, Nx = Nx, Dx = Dx, lx = lx, px = px,
            qx = qx, Lx = Lx, Tx = Tx, ex = ex)

# ======================== #
# PRODUCT LIMIT ESTIMATOR:
# ======================== #
n <- length(ageLast)

# create identities and age list:
if (is.null(ageFirst)) {
  allAges <- sort(ageLast)
  allAgesId <- departType
} else {
  idAgeFirst <- which(ageFirst > min(ageLast))
  allAges <- c(ageFirst[idAgeFirst], ageLast)
  allAgesId <- c(rep("F", length(idAgeFirst)), departType)
}

ageTypes <- unique(allAgesId)
ntypes <- length(ageTypes)
idsort <- sort.int(allAges, index.return = TRUE)$ix
allAges <- allAges[idsort]
allAgesId <- allAgesId[idsort]
nAllAges <- length(allAges)
unAllAges <- unique(allAges)
nages <- length(unAllAges)
ageNames <- as.character(unAllAges)

# Count by type:
recTab <- matrix(0, nages, ntypes, dimnames = list(ageNames, ageTypes))
for (at in ageTypes) {
  idtemp <- rep(0, nAllAges)
  idtemp[which(allAgesId == at)] <- 1
  ttemp <- table(allAges, idtemp)
  temp <- c(ttemp[, 2])
  recTab[, at] <- temp
}

# Cumulative tables:
cumTab <- apply(recTab, 2, function(xx) {
  rev(cumsum(rev(xx)))
})

Nx <- cumTab[, "D"] + cumTab[, "C"] - cumTab[, "F"]
Dx <- recTab[, "D"]

idDead <- which(recTab[, "D"] > 0)
ple <- cumprod(1 - c(Dx / Nx)[idDead])

# subset ages:
pleAges <- unAllAges[idDead]
# Add age 0:
if (pleAges[1] > 0) {
  pleAges <- c(0, pleAges)
  ple <- c(1, ple)
}

# Fill up table:
pleTab <- data.frame(Ages = pleAges, ple = ple)

# ============================ #
# OLD PRODUCT LIMIT ESTIMATOR:
# ============================ #
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
pleTab2 <- data.frame(Ages = agev, ple = ple)


# ====== #
# PLOTS: 
# ====== #
]