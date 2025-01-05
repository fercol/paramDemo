# ============================ CODE METADATA ================================= #
# PACKAGE: paramDemo
# AUTHOR: Fernando Colchero
# DATE CREATED: 2020-06-01
# DATE MODIFIED: 2023-11-18
# DESCRIPTION: Functions to extract demographic information from parametric
#              mortality and Fertility models, summary statistics (e.g.
#              ageing rates, life expectancy, lifespan equality, etc.), 
#              life table and product limit estimator calculation from 
#              census data.
# COMMENTS: a) Created a function to calculate CIs of life table elements and
#           plotting functions for life tables and life table CIs.
#           b) Created functions to calculate product limit 
#           estimators and a function to calculate their CIs using bootstrap.
# ============================ START CODE ==================================== #
# ============================= #
# ==== INTERNAL FUNCTIONS: ====
# ============================= #
# ---------------------------------------- #
# VERIFY AND EXTRACT PARAMETER ATTRIBUTES:
# ---------------------------------------- #
# Verify survival model and shape:
.VerifySurvMod <- function(model, shape) {
  if (!model %in% c("EX", "GO", "WE", "LO")) {
    stop("Wrong 'model' specification.\n",
         "       Alternatives are 'EX', 'GO', 'WE', or 'LO'.", call. = FALSE)
  }
  if (!shape %in% c("simple", "Makeham", "bathtub")) {
    stop("Wrong 'shape' specification.\n",
         "       Alternatives are 'simple', 'Makeham', or 'bathtub'.",
         call. = FALSE)
  }
  if (model == "EX" & shape != "simple") {
    warning("Only shape 'simple' can be used with model 'EX'", call. = FALSE)
  }
}

# Verify fertility model:
.VerifyFertMod <- function(modelFert) {
  fMods <- c("quadratic", "PeristeraKostaki", "ColcheroMuller", 
             "Hadwiger", "gamma", "beta", "skewNormal", 'gammaMixture',
             "HadwigerMixture", "skewSymmetric", "skewLogistic")
  if (!modelFert %in% fMods) {
    stop("Wrong 'modelFert' specification. See help file for fertility models.",
         call. = FALSE)
  }
}

# Survival:
.SetTheta  <- function(theta, model = "GO", shape = "simple") {
  if (is.null(theta)) {
    stop("Missing 'theta' parameter vector or matrix.\n", call. = FALSE)
  } 
  if (model == "EX") {
    nTh <- 1
    nameTh <- "b0"
    lowTh <- 0
  } else if (model == "GO") {
    nTh <- 2
    nameTh <- c("b0", "b1")
    lowTh <- c(-Inf, 0)
  } else if (model == "WE") {
    nTh <- 2
    nameTh <- c("b0", "b1")
    lowTh <- c(0, 0)
  } else if (model == "LO") {
    nTh <- 3
    nameTh <- c("b0", "b1", "b2")
    lowTh <- c(-Inf, 0, 0)
  }
  if (shape == "Makeham") {
    nTh <- nTh + 1
    nameTh <- c("c", nameTh)
    lowTh <- c(0, lowTh)
    if (theta[3] < 0 & model == "GO") {
      lowTh[3] <- -Inf
    }
  } else if (shape == "bathtub") {
    nTh <- nTh + 3
    nameTh <- c("a0", "a1", "c", nameTh)
    lowTh <- c(-Inf, 0, 0, lowTh)
  }
  if (is.matrix(theta)) {
    nthUser <- ncol(theta)
    clTh <- "matrix"
    stTh <- "columns"
  } else if (is.numeric(theta)) {
    nthUser <- length(theta)
    clTh <- "vector"
    stTh <- "elements"
  } else {
    stop("Parameters theta should either be of class matrix\nor a numeric vector.\n", call. = FALSE)
  }
  # Verify dimensions of theta vector or matrix:
  if (nthUser > nTh) {
    stop(sprintf("The theta %s should have %s %s.\n", clTh, nTh, stTh),
         call. = FALSE)
  } else {
    names(theta) <- nameTh
  }

  # Verify that all thetas are within the parameters support:
  if (is.matrix(theta)) {
    THELOW <- all(sapply(1:nTh, function(ti) {
      tl <- all(theta[, ti] >= lowTh[ti])
    }))
  } else {
    THELOW <- all(theta >= lowTh)
  }
  if(!THELOW) {
    stop(sprintf("Some theta parameters are below their lower bound.\n %s.\n",
                 paste(sprintf("min(%s) = %s", nameTh, lowTh), 
                       collapse = ", ")),
         call. = FALSE)
    
  }
  
  defaultTheta  <- list(theta = theta, p = nTh, name = nameTh,
                        low = lowTh)
  attr(defaultTheta, "model") = model
  attr(defaultTheta, "shape") = shape
  return(defaultTheta)
}

# Define beta:
.DefineBeta <- function(modelFert = "quadratic") {
  if (modelFert == "quadratic") {
    nBe <- 3
    lowBe <- rep(0, 3)
    uppBe <- rep(Inf, 3)
  } else if (modelFert == "PeristeraKostaki") {
    nBe <- 4
    lowBe <- c(0, 0, 0, 0)
    uppBe <- rep(Inf, 4)
  } else if (modelFert == "ColcheroMuller") {
    nBe <- 4
    lowBe <- c(0, 0, 0, -Inf)
    uppBe <- rep(Inf, 4)
  } else if (modelFert == "Hadwiger") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
    uppBe <- rep(Inf, 3)
  } else if (modelFert == "gamma") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
    uppBe <- rep(Inf, 3)
  } else if (modelFert == "beta") {
    nBe <- 5
    lowBe <- rep(0, 5)
    uppBe <- rep(Inf, 5)
  } else if (modelFert == "skewNormal") {
    nBe <- 4
    lowBe <- rep(0, 4)
    uppBe <- rep(Inf, 4)
  } else if (modelFert == "gammaMixture") {
    nBe <- 6
    lowBe <- rep(0, 6)
    uppBe <- c(Inf, 1, rep(Inf, 4))
  } else if (modelFert == "HadwigerMixture") {
    nBe <- 5
    lowBe <- rep(0, 5)
    uppBe <- c(1, rep(Inf, 4))
  } else if (modelFert == "skewSymmetric") {
    nBe <- 5
    lowBe <- c(rep(0, 4), -Inf)
    uppBe <- rep(Inf, 5)
  } else if (modelFert == "skewLogistic") {
    nBe <- 5
    lowBe <- c(rep(0, 4), -Inf)
    uppBe <- rep(Inf, 5)
  }
  nameBe <- sprintf("b%s", 1:nBe - 1)
  defaultBeta  <- list(p = nBe, name = nameBe, low = lowBe, upp = uppBe)
  return(defaultBeta)
}

# Fertility:
.SetBeta <- function(beta, modelFert = "quadratic") {
  if (is.null(beta)) {
    stop("Missing 'beta' parameter vector or matrix.\n", call. = FALSE)
  } 
  defBeta <- .DefineBeta(modelFert = modelFert)
  
  if (is.matrix(beta)) {
    nbeUser <- ncol(beta)
    clBe <- "matrix"
    stBe <- "columns"
  } else if (is.numeric(beta)) {
    nbeUser <- length(beta)
    clBe <- "vector"
    stBe <- "elements"
  } else {
    stop("The beta parameters should either be of class matrix", 
         " or a numeric vector.\n", call. = FALSE)
  }
  if (nbeUser != defBeta$p) {
    stop(sprintf("The beta %s should have %s %s.\n", clBe, defBeta$p, stBe),
         call. = FALSE)
  } else {
    if (is.matrix(beta)) {
      colnames(beta) <- defBeta$name
    } else {
      names(beta) <- defBeta$name
    }
  }
  # check if beta parameters conform to their support:
  if (is.matrix(beta)) {
    BETLOW <- all(sapply(1:defBeta$p, function(bi) {
      bl <- all(beta[, bi] >= defBeta$low[bi])
    }))
  } else {
    BETLOW <- all(beta >= defBeta$low)
  }
  if (is.matrix(beta)) {
    BETUPP <- all(sapply(1:defBeta$p, function(bi) {
      bl <- all(beta[, bi] <= defBeta$upp[bi])
    }))
  } else {
    BETUPP <- all(beta <= defBeta$upp)
  }
  if(!BETLOW) {
    stop(sprintf("Some beta parameters are below their lower bound.\n %s.\n",
                 paste(sprintf("min(%s) = %s", defBeta$name, defBeta$low), 
                       collapse = ", ")),
         call. = FALSE)
    
  }
  if(!BETUPP) {
    stop(sprintf("Some beta parameters are above their upper bound.\n %s.\n",
                 paste(sprintf("min(%s) = %s", defBeta$name, defBeta$low), 
                       collapse = ", ")),
         call. = FALSE)
    
  }
  defaultBeta  <- list(beta = beta, p = defBeta$p, name = defBeta$name,
                       low = defBeta$low, upp = defBeta$upp)
  attr(defaultBeta, "model") = modelFert
  return(defaultBeta)
}

# ---------------------- #
# DEMOGRAPHIC FUNCTIONS:
# ---------------------- #
# BASIC MORTALITY, CUMULATIVE HAZARDS AND Fertility FUNCTIONS:
# a) Mortality:
.DefineMortMatrix <- function(model = "GO", shape = "simple") {
  if (model == "EX") {
    mortfun <- function(theta, x) c(theta) * rep(1, length(x))
  } else if (model == "GO") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta[ ,"b0"] + theta[, "b1"] * x)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * x)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] +
          exp(theta[, "b0"] + theta[, "b1"] * x)
      }
    }
  } else if (model == "WE") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        theta[, "b0"] * theta[, "b1"]^theta[, "b0"] *
          (x + 0.003)^(theta[, "b0"] - 1)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + theta[, "b0"] * theta[, "b1"]^theta[, "b0"] *
          (x + 0.003)^(theta[, "b0"] - 1)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] +
          theta[, "b0"] * theta[, "b1"]^theta[, "b0"] *
          (x + 0.003)^(theta[, "b0"] - 1)
      }
    }
  } else if (model == "LO") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta[, "b0"] + theta[, "b1"] * x) /
          (1 + theta[, "b2"] * exp(theta[, "b0"]) /
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * x) /
          (1 + theta[, "b2"] * exp(theta[, "b0"]) /
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] +
          exp(theta[, "b0"] + theta[, "b1"] * x) /
          (1 + theta[, "b2"] * exp(theta[, "b0"]) /
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    }
  }
  return(mortfun)
}

.DefineMortNumeric <- function(model = "GO", shape = "simple") {
  if (model == "EX") {
    mortfun <- function(theta, x) c(theta) * rep(1, length(x))
  } else if (model == "GO") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta["b0"] + theta["b1"] * x)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + exp(theta["b0"] + theta["b1"] * x)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
          exp(theta["b0"] + theta["b1"] * x)
      }
    }
  } else if (model == "WE") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        theta["b0"] * theta["b1"]^theta["b0"] *
          (x + 0.003)^(theta["b0"] - 1)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + theta["b0"] * theta["b1"]^theta["b0"] *
          (x + 0.003)^(theta["b0"] - 1)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
          theta["b0"] * theta["b1"]^theta["b0"] *
          (x + 0.003)^(theta["b0"] - 1)
      }
    }
  } else if (model == "LO") {
    if (shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta["b0"] + theta["b1"] * x) /
          (1 + theta["b2"] * exp(theta["b0"]) /
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + exp(theta["b0"] + theta["b1"] * x) /
          (1 + theta["b2"] * exp(theta["b0"]) /
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
          exp(theta["b0"] + theta["b1"] * x) /
          (1 + theta["b2"] * exp(theta["b0"]) /
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    }
  }
  return(mortfun)
}

# b) Cummulative hazard:
.DefineCumHazMatrix <- function(model = "GO", shape = "simple") {
  if (model == "EX") {
    cumhazfun <- function(theta, x) c(theta) * x
  } else if (model == "GO") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        exp(theta[, "b0"]) / theta[, "b1"] *
          (exp(theta[, "b1"] * x) - 1)
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] *
          (exp(theta[, "b1"] * x) - 1)
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
          theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] *
          (exp(theta[, "b1"] * x) - 1)
      }
    }
  } else if (model == "WE") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        (theta[, "b1"] * x)^theta[, "b0"]
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + (theta[, "b1"] * x)^theta[, "b0"]
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
          theta[, "c"] * x + (theta[, "b1"] * x)^theta[, "b0"]
      }
    }
  } else if (model == "LO") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        log(1 + theta[, "b2"] * exp(theta[, "b0"]) / theta[, "b1"] *
              (exp(theta[, "b1"] * x) - 1)) * (1 / theta[, "b2"])
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + log(1 + theta[, "b2"] * exp(theta[, "b0"]) /
                                 theta[, "b1"] *
                                 (exp(theta[, "b1"] * x) - 1)) *
          (1 / theta[, "b2"])
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
          theta[, "c"] * x + log(1 + theta[, "b2"] *
                                   exp(theta[, "b0"]) / theta[, "b1"] *
                                   (exp(theta[, "b1"] * x) - 1)) *
          (1 / theta[, "b2"])
      }
    }
  }
  return(cumhazfun)
}

.DefineCumHazNumeric <- function(model = "GO", shape = "simple") {
  if (model == "EX") {
    cumhazfun <- function(theta, x) c(theta) * x
  } else if (model == "GO") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        exp(theta["b0"]) / theta["b1"] *
          (exp(theta["b1"] * x) - 1)
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + exp(theta["b0"]) / theta["b1"] *
          (exp(theta["b1"] * x) - 1)
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
          theta["c"] * x + exp(theta["b0"]) / theta["b1"] *
          (exp(theta["b1"] * x) - 1)
      }
    }
  } else if (model == "WE") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        (theta["b1"] * x)^theta["b0"]
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + (theta["b1"] * x)^theta["b0"]
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
          theta["c"] * x + (theta["b1"] * x)^theta["b0"]
      }
    }
  } else if (model == "LO") {
    if (shape == "simple") {
      cumhazfun <- function(theta, x) {
        log(1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
              (exp(theta["b1"] * x) - 1)) * (1 / theta["b2"])
      }
    } else if (shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + log(1 + theta["b2"] * exp(theta["b0"]) /
                               theta["b1"] *
                               (exp(theta["b1"] * x) - 1)) *
          (1 / theta["b2"])
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
          theta["c"] * x + log(1 + theta["b2"] *
                                 exp(theta["b0"]) / theta["b1"] *
                                 (exp(theta["b1"] * x) - 1)) *
          (1 / theta["b2"])
      }
    }
  }
  return(cumhazfun)
}

# c) Fertility:
.DefineFertilityNumeric <- function(modelFert = "quadratic") {
  if (modelFert == "quadratic") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2)
      return(fert)
    }
  } else if (modelFert == "PeristeraKostaki") {
    fertfun <- function(beta, x) {
      be1 <- x * 0 + beta["b1"]
      be1[which(x > beta["b3"])] <- beta["b2"]
      fert <- beta["b0"] * exp(-((x - beta["b3"]) / be1)^2)
      return(fert)
    }
  } else if (modelFert == "ColcheroMuller") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2 +
                                 beta["b3"] * 1/(x + 1))
      return(fert)
    }
  } else if (modelFert == "Hadwiger") {
    fertfun <- function(beta, x) {
      fert <- (beta["b0"] * beta["b1"])/beta["b2"] * 
        (beta["b2"]/x)^(3/2) * exp(-beta["b1"]^2 * (beta["b2"] / x + 
                                                      x / beta["b2"] - 2))
      return(fert)
    }
  } else if (modelFert == "gamma") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * dgamma(x, shape = beta["b1"], rate = beta["b2"])
      return(fert)
    }
  } else if (modelFert == "beta") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * ((x - beta["b3"])^(beta["b1"] - 1) * 
                              (beta["b4"] - x)^(beta["b2"] - 1)) / 
        ((beta["b4"] - beta["b3"])^(beta["b1"] + beta["b2"] - 1) * 
           beta(beta["b1"], beta["b2"]))
      return(fert)
    }
  } else if (modelFert == "skewNormal") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1/beta["b1"] * 
        dnorm((x - beta["b2"]) / beta["b1"]) * 
        pnorm(beta["b3"] * ((x - beta["b2"]) / beta["b1"]))
      return(fert)
    }
  } else if (modelFert == "gammaMixture") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * (beta["b1"] * 
        dgamma(x, shape = beta["b2"], rate = beta["b3"]) +
      (1 - beta["b1"]) * dgamma(x, shape = beta["b4"], rate = beta["b5"]))
      return(fert)
    }
  } else if (modelFert == "HadwigerMixture") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * (beta["b1"]/beta["b2"] * 
        (beta["b2"]/x)^(3/2) * exp(-beta["b1"]^2 * (beta["b2"] / x + 
                                                      x / beta["b2"] - 2))) +
        (1 - beta["b0"]) * (beta["b3"]/beta["b4"] * 
                              (beta["b4"]/x)^(3/2) * 
                              exp(-beta["b3"]^2 * (beta["b4"] / x + 
                                                     x / beta["b4"] - 2)))
      return(fert)
    }
  } else if (modelFert == "skewSymmetric") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1/beta["b1"] * 
        dnorm((x - beta["b2"]) / beta["b1"]) * 
        pnorm(beta["b3"] * ((x - beta["b2"]) / beta["b1"]) + 
                beta["b4"] * ((x - beta["b2"]) / beta["b1"])^3)
      return(fert)
    }
  } else if (modelFert == "skewLogistic") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * 2 * 1 / beta["b1"] * 
        ((exp(-(x - beta["b2"]) / beta["b1"])) / 
           ((1 + exp(-(x - beta["b2"]) / beta["b1"]))^2 *
              (1 + exp(-beta["b3"] * (x - beta["b2"]) / beta["b1"] - 
                         beta["b4"] * ((x - beta["b2"]) / beta["b1"])^3))))
      return(fert)
    }
  }
  return(fertfun)
}

.DefineFertilityMatrix <- function(modelFert = "quadratic") {
  if (modelFert == "quadratic") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2)
      return(fert)
    }
  } else if (modelFert == "PeristeraKostaki") {
    fertfun <- function(beta, x) {
      be1 <- x * 0 + beta[, "b1a"]
      be1[which(x > beta[, "b2"])] <- beta[, "b1b"]
      fert <- beta[, "b0"] * exp(-((x - beta[, "b2"]) / be1)^2)
      return(fert)
    }
  } else if (modelFert == "ColcheroMuller") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2 +
                                   beta[, "b3"] * 1/(x + 1))
      return(fert)
    }
  } else if (modelFert == "Hadwiger") {
    fertfun <- function(beta, x) {
      fert <- (beta[, "b0"] * beta[, "b1"])/beta[, "b2"] * 
        (beta[, "b2"]/x)^(3/2) * exp(-beta[, "b1"]^2 * 
                                       (beta[, "b2"] / x + 
                                          x / beta[, "b2"] - 2))
      return(fert)
    }
  } else if (modelFert == "gamma") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * dgamma(x, shape = beta[, "b1"], 
                                    rate = beta[, "b2"])
      return(fert)
    }
  } else if (modelFert == "beta") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * ((x - beta[, "b3"])^(beta[, "b1"] - 1) * 
                                (beta[, "b4"] - x)^(beta[, "b2"] - 1)) / 
        ((beta[, "b4"] - beta[, "b3"])^(beta[, "b1"] + beta[, "b2"] - 1) * 
           beta(beta[, "b1"], beta[, "b2"]))
      return(fert)
    }
  } else if (modelFert == "skewNormal") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1/beta[, "b1"] * 
        dnorm((x - beta[, "b2"]) / beta[, "b1"]) * 
        pnorm(beta[, "b3"] * ((x - beta[, "b2"]) / beta[, "b1"]))
      return(fert)
    }
  } else if (modelFert == "gammaMixture") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * (beta[, "b1"] * 
                                dgamma(x, shape = beta[, "b2"], 
                                       rate = beta[, "b3"]) +
                                (1 - beta[, "b1"]) * 
                                dgamma(x, shape = beta[, "b4"], 
                                       rate = beta[, "b5"]))
      return(fert)
    }
  } else if (modelFert == "HadwigerMixture") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * (beta[, "b1"]/beta[, "b2"] * 
                                (beta[, "b2"]/x)^(3/2) * 
                                exp(-beta[, "b1"]^2 * 
                                      (beta[, "b2"] / x + 
                                         x / beta[, "b2"] - 2))) +
        (1 - beta[, "b0"]) * (beta[, "b3"]/beta[, "b4"] * 
                                (beta[, "b4"]/x)^(3/2) * 
                                exp(-beta[, "b3"]^2 * (beta[, "b4"] / x + 
                                                         x / beta[, "b4"] - 2)))
      return(fert)
    }
  } else if (modelFert == "skewSymmetric") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1/beta[, "b1"] * 
        dnorm((x - beta[, "b2"]) / beta[, "b1"]) * 
        pnorm(beta[, "b3"] * ((x - beta[, "b2"]) / beta[, "b1"]) + 
                beta[, "b4"] * ((x - beta[, "b2"]) / beta[, "b1"])^3)
      return(fert)
    }
    
  } else if (modelFert == "skewLogistic") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * 2 * 1 / beta[, "b1"] * 
        ((exp(-(x - beta[, "b2"]) / beta[, "b1"])) / 
           ((1 + exp(-(x - beta[, "b2"]) / beta[, "b1"]))^2 *
              (1 + exp(-beta[, "b3"] * (x - beta[, "b2"]) / beta[, "b1"] - 
                         beta[, "b4"] * ((x - beta[, "b2"]) / 
                                           beta[, "b1"])^3))))
      return(fert)
    }
    
  }
  return(fertfun)
}

# ------------------------ #
# NON-PARAMETRIC SURVIVAL:
# ------------------------ #
# Life table:
.CalcLT <- function(ageLast, ageFirst, departType, dx) {
  # Unit age vector for that sex:
  agev <- seq(from = 0, to = ceiling(max(ageLast[which(departType == "D")])), 
              by = dx)
  nage <- length(agev)
  
  # Outputs:
  Nx <- Dx <- ax <- rep(0, nage)
  for (ix in 1:nage) {
    # A) EXPOSURES:
    # Find how many entered the interval (including truncated):
    idNx <- which(ageFirst < agev[ix] + dx & ageLast >= agev[ix])
    nNx <- length(idNx)
    
    # Extract ages and departType:
    xFirst <- ageFirst[idNx]
    xLast <- ageLast[idNx]
    dType <- departType[idNx]
    
    # Index for individuals dying within interval:
    idDx <- which(xLast < agev[ix] + dx & dType == "D")
    nDx <- length(idDx)
    
    # Index of truncated in interval:
    idtr <- which(xFirst >= agev[ix])
    
    # Index of censored in the interval:
    idce <- which(xLast < agev[ix] + dx & dType == "C")
    
    # Porportion lived within interval:
    intr <- rep(0, nNx)
    ince <- rep(dx, nNx)
    intr[idtr] <- xFirst[idtr] - agev[ix]
    ince[idce] <- xLast[idce] - agev[ix]
    lived <- (ince - intr) / dx
    
    # Fill in Nx:
    Nx[ix] <- sum(lived)
    
    # B) DEATHS:
    # Fill in Dx:
    Dx[ix] <- nDx
    
    
    # C) PROPORTION LIVED BY THOSE THAT DIED IN INTERVAL:
    if (Dx[ix] > 1) {
      ylived <- xLast[idDx] - agev[ix]
      ax[ix] <- sum(ylived) / dx / nDx
    } else {
      ax[ix] <- 0
    }
  }
  
  # Age-specific mortality probability:
  qx <- Dx / Nx
  qx[which(is.na(qx))] <- 0
  
  # Age-specific survival probability:
  px <- 1 - qx
  
  # Survivorship (or cumulative survival):
  lx <- c(1, cumprod(px))[1:nage]
  
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
  
  # Life-table:
  lt <- cbind(Ages = agev, Nx = Nx, Dx = Dx, lx = lx, px = px,
              qx = qx, Lx = Lx, Tx = Tx, ex = ex)
  return(lt)
}

# Product limit estimator:
.CalcPLE <- function(ageLast, ageFirst, departType) {
  # Sort ages:
  idsort <- sort.int(ageLast, index.return = TRUE)$ix
  
  # Create new age vector:
  agev <- unique(ageLast[idsort])
  
  # Number of unique ages:
  nage <- length(agev)
  
  # Cx and delta x vector:
  Cx <- rep(0, nage)
  delx <- rep(0, nage)
  
  # Fill up Cx and delta:
  for (ii in 1:nage) {
    idNx <- which(ageFirst <= agev[ii] & ageLast >= agev[ii])
    Cx[ii] <- length(idNx) / nage
    idd <- which(ageLast == agev[ii] & departType == "D")
    delx[ii] <- length(idd)
  }
  
  # Calculate product limit estimator:
  ple <- cumprod((1 - 1 / (nage * Cx))^delx)
  
  # Add age 0:
  if (agev[1] > 0) {
    agev <- c(0, agev)
    ple <- c(1, ple)
  }
  # Create data frame:
  pleTab <- data.frame(Ages = agev, ple = ple)
  
  return(pleTab)
}

# ------------- #
# AGEING RATES:
# ------------- #
# Mortality (actuarial) ageing rates:
.DefineARmort <- function(model = "GO", shape = "simple") {
  # -------------------------------------
  # Exponential (i.e. constat mortality):
  # -------------------------------------
  if (model == "EX") {
    # ------------ #
    # Exponential:
    # ------------ #
    ageingRate <- function(theta, x) rep(0, length(x))
    # --------- #
    # Gompertz:
    # --------- #
  } else if (model == "GO") {
    if (shape == "simple") {
      ageingRate <- function(theta, x) rep(theta["b1"], length(x))
    } else if (shape == "Makeham") {
      ageingRate <- function(theta, x) {
        theta["b1"] * exp(theta["b0"] + theta["b1"] * x) /
          (theta["c"] + exp(theta["b0"] + theta["b1"] * x))
      }
    } else {
      ageingRate <- function(theta, x) {
        (-theta["a1"] * exp(theta["a0"] - theta["a1"] * x) +
           theta["b1"] * exp(theta["b0"] + theta["b1"] * x)) /
          (exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
             exp(theta["b0"] + theta["b1"] * x))
      }
    }
    # -------- #
    # Weibull:
    # -------- #
  } else if (model == "WE") {
    if (shape == "simple") {
      ageingRate <- function(theta, x) {
        (theta["b1"] * (theta["b0"] - 1)) / (theta["b1"] * x)
      }
    } else if (shape == "Makeham") {
      ageingRate <- function(theta, x) {
        (theta["b0"] * theta["b1"]^(theta["b0"]) * (theta["b0"] - 1) *
           x^(theta["b0"] - 2)) /
          (theta["c"] + theta["b0"] * theta["b1"]^(theta["b0"]) *
             x^(theta["b0"] - 1))
      }
    } else {
      ageingRate <- function(theta, x) {
        (-theta["a1"] * exp(theta["a0"] - theta["a1"] * x) +
           theta["b0"] * theta["b1"]^(theta["b0"]) * (theta["b0"] - 1) *
           x^(theta["b0"] - 2)) /
          (exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
             theta["b0"] * theta["b1"]^(theta["b0"]) * x^(theta["b0"] - 1))
      }
    }
    # --------- #
    # Logistic:
    # --------- #
  } else {
    if (shape == "simple") {
      ageingRate <- function(theta, x) {
        theta["b1"] - theta["b2"] * exp(theta["b0"] + theta["b1"] * x) /
          (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
             (exp(theta["b1"] * x) - 1))
      }
    } else if (shape == "Makeham") {
      ageingRate <- function(theta, x) {
        exp(theta["b0"] + theta["b1"] * x) *
          (theta["b1"] - theta["b2"] * exp(theta["b0"])) /
          ((1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
              (exp(theta["b1"] * x) - 1))^2 *
             (theta["c"] + exp(theta["b0"] + theta["b1"] * x) /
                (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
                   (exp(theta["b1"] * x) - 1))))
      }
    } else {
      ageingRate <- function(theta, x) {
        (-theta["a1"] * exp(theta["a0"] - theta["a1"] * x) +
           exp(theta["b0"] + theta["b1"] * x) *
           (theta["b1"] * (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
                             (exp(theta["b1"] * x) - 1)) -
              theta["b2"] * exp(theta["b0"] + theta["b1"] * x)) /
           (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
              (exp(theta["b1"] * x) - 1))^2) * 1 /
          (exp(theta["a0"] - theta["a1"] * x) +
             theta["c"] + exp(theta["b0"] + theta["b1"] * x) /
             (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
                (exp(theta["b1"] * x) - 1)))
      }
    }
  }
  return(ageingRate)
}

# ------------------------------------ #
# MORTALITY RATE DOUBLING TIME (MRDT):
# ------------------------------------ #
.CalcMRDT <- function(theta, x, model = "GO", shape = "simple",
                      checkTheta = TRUE, error = 1e-10) {

  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)

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
      MRDT <- log((2 - 2 * g) / (1 - g * (exp(theta["b1"] * x) + 1))) /
        theta["b1"]
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
      MRDT <- log((M * (1 - g)) / (1 - M * g * exp(theta["b1"] * x))) /
        theta["b1"]
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

# ----------------------------- #
# SUMMARY STATISTICS FUNCTIONS:
# ----------------------------- #
# Life expectancy, lifespan inequality, lifespan equality, Gini, coeff. of
# variation:
.CalcSummStatsMort <- function(theta, x, dx, survFunc) {
  Sx <- survFunc(theta, x)
  Sx <- Sx / Sx[1]
  idd <- which(Sx > 0)
  Sx <- Sx[idd]
  dS <- -diff(Sx)
  dS <- dS / sum(dS)
  nx <- length(Sx)
  Asx <- sapply(c(0.5, 0.2, 0.05), function(sxb) {
    .FindMaxAge(theta, bound = sxb, survFunc = survFunc)
  })
  names(Asx) <- sprintf("AgeSx=%s", c(0.5, 0.2, 0.05))
  Ex <- sum(Sx * dx) / Sx[1]
  Hx <- -sum(((Sx * log(Sx))[-1] + (Sx * log(Sx))[-nx]) / 2 * dx) / Ex
  Kx <- -log(Hx)
  Gx <- 1 - 1 / sum((Sx[-1] + Sx[-nx]) / 2 * dx) *
    sum((Sx[-1]^2 + Sx[-nx]^2) / 2 * dx)
  x <- (x - x[1])[idd]
  CV <- sqrt(sum((x[-length(x)] + dx/2 - Ex)^2 * dS)) / Ex

  return(c(Asx, lifeExp = Ex, lifespIneq = Hx, lifespEqual = Kx,
           Gini = Gx, CoeffVar = CV))
}

# ------------------------------------- #
# FIND UPPER BOUND AGE FOR INTEGRATION:
# ------------------------------------- #
.FindMaxAge <- function(theta, bound = 0.0001, error = 10^(-5), survFunc) {
  xTest <- c(0, 10^(0:10))
  SxTest <- survFunc(theta, xTest)
  idBound <- which(SxTest < bound)[1]
  SxDiff <- abs(SxTest[idBound - 1] - bound)
  icount <- 0
  while(SxDiff > error) {
    icount <- icount + 1
    xTest <- seq(xTest[idBound - 1], xTest[idBound], length = 100)
    SxTest <- survFunc(theta, xTest)
    idBound <- which(SxTest < bound)[1]
    SxDiff <- abs(SxTest[idBound - 1] - bound)
  }
  return(ceiling(xTest[idBound - 1]))
}

# ---------------------------- #
# TRANSITION MATRIX FUNCTIONS:
# ---------------------------- #
# Function to construct Leslie matrix:
.FillLeslieMatr <- function(p, f, n) {
  idcol <- 1:n
  idrow <- c(2:n, n)
  idpa <- (idcol - 1) * n + idrow
  Aa <- matrix(0, n, n)
  Aa[1, ] <- f
  Aa[idpa] <- p
  return(Aa)
}

# ========================= #
# ==== USER FUNCTIONS: ====
# ========================= #
# --------------------------- #
# BASIC PARAMETRIC FUNCTIONS: 
# --------------------------- #
# Calculate mortality:
CalcMort <- function(theta, x, model = "GO", shape = "simple",
                     checkTheta = TRUE) {

  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)

  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }

  # a) Mortality function:
  .CalcMort <- function(theta, ...) UseMethod(".CalcMort")
  .CalcMort.matrix <- .DefineMortMatrix(model = model, shape = shape)
  .CalcMort.numeric <- .DefineMortNumeric(model = model, shape = shape)

  # b) Calculate mortality:
  mort <- .CalcMort(theta, x)

  return(mort)
}

# Calculate Survival:
CalcSurv <- function(theta, x, model = "GO", shape = "simple",
                     checkTheta = TRUE) {

  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)

  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }

  # a) Cummulative hazard:
  .CalcCumHaz <- function(theta, ...) UseMethod(".CalcCumHaz")
  .CalcCumHaz.matrix <- .DefineCumHazMatrix(model = model, shape = shape)
  .CalcCumHaz.numeric <- .DefineCumHazNumeric(model = model, shape = shape)

  # b) Survival:
  .CalcSurv <- function(theta, x) {
    exp(-.CalcCumHaz(theta = theta, x = x))
  }

  # b) Calculate mortality:
  surv <- .CalcSurv(theta, x)

  return(surv)
}

# Calculate Fertility:
CalcFert <- function(beta, x, modelFert = "quadratic", checkBeta = TRUE) {

  # Verify fertility model:
  .VerifyFertMod(modelFert = modelFert)
  
  # Extract beta attributes:
  if (checkBeta) {
    betaAttr <- .SetBeta(beta, modelFert = modelFert)
    beta <- betaAttr$beta
  }

  # a) Fertility method:
  .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
  .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
  .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)

  # b) Calculate Fertility:
  fertfun <- .CalcFert(beta, x)
  return(fertfun)
}

# --------------------------------------------- #
# ---- GENERAL PARAMETRIC DEMOGRAPHIC FUNCTION: 
# --------------------------------------------- #
# Main demographic function:
CalcDemo <- function(theta = NULL, beta = NULL, x = NULL, dx = NULL, 
                     model = "GO", shape = "simple", modelFert = "quadratic", 
                     type = "both", minSx = 0.01, summarStats = TRUE, 
                     ageMatur = 0, maxAge = NULL, agesAR = NULL, 
                     SxValsAR = NULL) {
  if (!type %in% c("survival", "fertility", "both")) {
    stop("Wrong 'type' argument, values should be:
         'survival', 'fertility', or 'both'.",
         call. = FALSE)
  }
  # logical for survival:
  if (type %in% c("survival", "both")) {
    SURV <- TRUE
  } else {
    SURV <- FALSE
  }
  
  # Logical for fertility:
  if (type %in% c("fertility", "both")) {
    FERT <- TRUE
  } else {
    FERT <- FALSE
  }
  
  # =========== #
  # I) SURVIVAL: 
  # =========== #
  if (SURV) {
    # ------------------------------------------------ #
    # A) GENERAL SETUP AND VERIFICATION OF PARAMETERS:
    # ------------------------------------------------ #
    # Verify model and shape:
    .VerifySurvMod(model = model, shape = shape)
    
    # Extract theta attributes:
    thetaAttr <- .SetTheta(theta = theta, model = model, shape = shape)
    theta <- thetaAttr$theta
    
    # --------------------- #
    # B) PREPARE FUNCTIONS:
    # --------------------- #
    # a) Mortality function:
    .CalcMort <- function(theta, ...) UseMethod(".CalcMort")
    .CalcMort.matrix <- .DefineMortMatrix(model = model, shape = shape)
    .CalcMort.numeric <- .DefineMortNumeric(model = model, shape = shape)
    
    # b) Cummulative hazard:
    .CalcCumHaz <- function(theta, ...) UseMethod(".CalcCumHaz")
    .CalcCumHaz.matrix <- .DefineCumHazMatrix(model = model, shape = shape)
    .CalcCumHaz.numeric <- .DefineCumHazNumeric(model = model, shape = shape)
    
    # c) Survival:
    .CalcSurv <- function(theta, x) {
      exp(-.CalcCumHaz(theta = theta, x = x))
    }
    
    # d) Density:
    .CalcDens <- function(theta, x) {
      .CalcSurv(theta = theta, x = x) * .CalcMort(theta = theta, x = x)
    }
    
    # -------------------------------- #
    # C) CALCULATE DEMOGRAPHIC VALUES:
    # -------------------------------- #
    # Find upper age bound for plotting:
    if (is.null(x)) {
      if (is.null(maxAge)) {
        maxAge <- .FindMaxAge(theta, bound = minSx, survFunc = .CalcSurv)
      }
      if (is.null(dx)) dx <- 0.01
      x <- seq(0, maxAge, dx)
    }
    
    # mortality:
    mort <- .CalcMort(theta, x)
    
    # Cumulative hazards:
    cumhaz <- .CalcCumHaz(theta, x)
    
    # Survival:
    surv <- .CalcSurv(theta, x)
    
    # PDF:
    dens <- .CalcDens(theta, x)
    
    if (summarStats) {
      # Find upper age bound for integration:
      xMax <- .FindMaxAge(theta, survFunc = .CalcSurv)
      dxInt <- xMax / 5000
      xInt <- seq(0, xMax, dxInt)
      
      # Ageing rates:
      if (is.null(agesAR)) {
        if (is.null(SxValsAR)) SxValsAR <- c(0.5, 0.2, 0.05)
        Sx <- .CalcSurv(theta, xInt)
        agesAR <- sapply(SxValsAR, function(sxt) {
          idx <- which(abs(Sx - sxt) == min(abs(Sx - sxt)))
          return(xInt[idx])
        })
      } else {
        SxValsAR <- .CalcSurv(theta, agesAR)
      }
      ageingRates <- cbind(CalcAgeingRateMort(theta = theta, x = agesAR, 
                                          model = model, shape = shape), 
                           Surv = round(SxValsAR, 4))
      
      # Summary statistics for mortality:
      summStatsMort <- .CalcSummStatsMort(theta, xInt, dxInt, 
                                          survFunc = .CalcSurv)
      
      # Logical for calculation of summary statistics:
      calc <- TRUE
    } else {
      ageingRates <- NA
      summStatsMort <- NA
      calc <- FALSE
    }
    
    # List of outputs:
    survList <- list(functs = data.frame(age = x, mort = mort, surv = surv, 
                                         pdf = dens, cumhaz = cumhaz), 
                     summStats = list(calculated = calc, 
                                      ageingRates = ageingRates,
                                      summStatsMort = summStatsMort), 
                     settings = list(theta = theta, model = model,
                                     shape = shape), analyzed = TRUE)
  } else {
    survList <- list(analyzed = FALSE)
  } 
  
  # ============== #
  # II) FERTILITY: 
  # ============= #
  if (FERT) {
    # ------------------------------------------------ #
    # A) GENERAL SETUP AND VERIFICATION OF PARAMETERS:
    # ------------------------------------------------ #
    # Verify fertility model:
    .VerifyFertMod(modelFert = modelFert)
    
    # Extract beta attributes:
    betaAttr <- .SetBeta(beta = beta, modelFert = modelFert)
    beta <- betaAttr$beta
    
    # --------------------- #
    # B) PREPARE FUNCTIONS:
    # --------------------- #
    # Fertility function:
    .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
    .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)
    .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
    
    # -------------------------------- #
    # C) CALCULATE DEMOGRAPHIC VALUES:
    # -------------------------------- #
    if (is.null(x)) {
      if (is.null(maxAge)) {
        stop("Missing argument 'x' for vector of ages, or argument for maximum",
        " age 'maxAge'. Provide either.", call. = FALSE)
      } else {
        if (is.null(dx)) dx <- 0.01
        x <- seq(ageMatur, maxAge, dx)
      }
    }
    
    # Calculate age-specific fertility:
    fert <- .CalcFert(beta, x[which(x >= ageMatur)] - ageMatur) * dx
    
    # Summary 
    if (summarStats) {
      # Age at maximum fertility:
      ageMaxFert <- CalcAgeMaxFert(beta = beta, modelFert = modelFert, 
                                   ageMatur = ageMatur, maxAge = maxAge)
      calc <- TRUE
      
      # Ageing rates:
      if (is.null(agesAR)) {
        if (is.null(SxValsAR)) SxValsAR <- c(0.5, 0.2, 0.05)
        Sx <- .CalcSurv(theta, xInt)
        agesAR <- sapply(SxValsAR, function(sxt) {
          idx <- which(abs(Sx - sxt) == min(abs(Sx - sxt)))
          return(xInt[idx])
        })
      } else {
        SxValsAR <- .CalcSurv(theta, agesAR)
      }
      ageingRatesFert <- cbind(CalcAgeingRateFert(beta = beta, x = agesAR, 
                                              modelFert = modelFert, 
                                              ageMatur = ageMatur), 
                           Surv = round(SxValsAR, 4))
      
    } else {
      ageMaxFert <- NA
      ageingRatesFert <- NA
      calc <- FALSE
    }
    fertList <- list(functs = data.frame(age = x[which(x >= ageMatur)], 
                                         fert = fert),
                     summStats = list(calculated = calc, 
                                      ageingRates = ageingRatesFert,
                                      summStatsFert = ageMaxFert),
                     settings = list(beta = beta, modelFert = modelFert), 
                     analyzed = TRUE)
  } else {
    fertList <- list(analyzed = FALSE)
  }
  
  # ========================== #
  # III) CREATE OUTPUT OBJECT: 
  # ========================== #
  outList <- list(surv = survList, fert = fertList)
  class(outList) <- c("paramDemo", ifelse(SURV & FERT, "demoBoth", 
                                          ifelse(SURV, "demoSurv", "demoFert")))
  return(outList)
}

# Plotting demographic object:
plot.paramDemo <- function(x, demofun = "all", ...) {
  # User par settings:
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Additional arguments:
  args <- list(...)
  
  # Labels for demographic rates:
  if (inherits(x, "demoBoth")) {
    demoFun <- c("mort", "surv", "pdf", "fert")
  } else if (inherits(x, "demoSurv")) {
    demoFun <- c("mort", "surv", "pdf")
  } else {
    demoFun <- "fert"
    demofun <- "fert"
  }
  
  # number of demoFun:
  nDemo <- length(demoFun)
  
  # Check if demofun is properly provided:
  if (inherits(x, "demoBoth") | inherits(x, "demoSurv")) {
    if (!demofun %in% c(demoFun, "all")) {
      stop(sprintf("Argument 'demofun' incorrect, values should be\n%s", 
                   paste(sprintf("'%s'", c(demoFun, "all")), collapse = ", ")), 
           call. = FALSE)
    }
  } else {
    if (demofun != "fert") {
      warning("Argument 'demofun' incorrect, value should be 'fert'.", 
              call. = FALSE)
    }
  }
  
  # Plot names for demographic rates:
  demoFunNames <- c(mort = expression(paste("Mortality, ", mu, "(", 
                                            italic(x), ")")),
                    surv = expression(paste("Survival, ", italic(S), "(", 
                                            italic(x), ")")), 
                    pdf = expression(paste("PDF of ages at death, ", italic(f), 
                                           "(", italic(x), ")")),
                    fert = expression(paste("Fertility, ", italic(m), "(", 
                                            italic(x), ")")))
  # mar values:
  if ("mar" %in% names(args)) {
    mar <- args$mar
  } else {
    mar <- c(2, 4, 1, 1)
  }
  
  # Axis labels:
  if ("xlab" %in% names(args)) {
    xlab <- args$xlab
  } else {
    xlab <- expression(paste("Age, ", italic(x)))
  }
  if ("ylab" %in% names(args)) {
    ylab <- args$ylab
    if (demofun == "all" & length(ylab) != nDemo) {
      stop(sprintf("Argument 'ylab' should be of length %s for demofun = 'all'.\nNote functions plotted will be %s.", nDemo, 
                   paste(sprintf("'%s'", demoFun), collapse = ", ")), 
           call. = FALSE)
    }
  } else {
    if (demofun == "all") {
      ylab <- demoFunNames
    } else {
      ylab <- demoFunNames[demofun]
    }
  }
  
  # Type of line:
  if ("type" %in% names(args)) {
    type <- args$type
  } else {
    type <- "l"
  }
  
  # Color
  if ("col" %in% names(args)) {
    col <- args$col
  } else {
    col <- "#800026"
  }
  
  # Line width:
  if ("lwd" %in% names(args)) {
    lwd <- args$lwd
  } else {
    lwd <- 1
  }
  
  # Cex.lab:
  if ("cex.lab" %in% names(args)) {
    cex.lab <- args$cex.lab
  } else {
    cex.lab <- 1
  }
  
  # Cex.axis:
  if ("cex.axis" %in% names(args)) {
    cex.axis <- args$cex.axis
  } else {
    cex.axis <- 1
  }
  
  # Produce plots:
  if (demofun == "all") {
    # -------------- #
    # MULTIPLE PLOTS:
    # -------------- #
    # Plot layout settings:
    if (inherits(x, "demoBoth")) {
      # Layout matrix:
      laymat <- cbind(c(2, 4, 1), c(3, 5, 1))
      
      # Widths and heights:
      widths <- c(1, 1)
      heights <- c(1, 1, 0.15)
    } else if (inherits(x, "demoSurv")) {
      # Layout matrix:
      laymat <- matrix(c(2, 3, 4, 1), ncol = 1)
      
      # Widths and heights:
      widths <- 1
      heights <- c(1, 1, 1, 0.15)
    }
    
    
    # produce plot:
    layout(mat = laymat, widths = widths, heights = heights)
    
    # x-axis label:
    par(mar = mar * c(0, 1, 0, 1))
    plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    text(0.5, 0.5, xlab, cex = 1.5 * cex.lab)
    
    par(mar = mar) 
    for (dr in demoFun) {
      if (dr == "fert") {
        agev <- x$fert$functs$age
        ydr <- x$fert$functs$fert
      } else {
        agev <- x$surv$functs$age
        ydr <- x$surv$functs[[dr]]
      }
      nage <- length(agev)
      if ("xlim" %in% names(args)) {
        xlim <- args$xlim
      } else {
        xlim <- range(agev) * c(0, 1.1)
      }
      if ("ylim" %in% names(args)) {
        ylim <- args$ylim
      } else {
        if (dr == "surv") {
          ylim <- c(0, 1)
        } else {
          ylim <- c(0, max(ydr))
        }
      }
      plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
      lines(agev, ydr, type = type, col = col, lwd = lwd)
      Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
      Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
      mtext(ylab[dr], side = 2, line = mar[2] / 2, cex = cex.lab)
    }    
  } else {
    # ------------ #
    # SINGLE PLOT:
    # ------------ #
    dr <- demofun
    # Layout matrix:
    laymat <- matrix(c(2, 1), nrow = 2, ncol = 1)
    
    # Widths and heights:
    widths <- c(1)
    heights <- c(1, 0.15)
    
    # produce plot:
    layout(mat = laymat, widths = widths, heights = heights)
    
    # x-axis label:
    par(mar = mar * c(0, 1, 0, 1))
    plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    # text(0.5, 0.5, xlab, cex = 1.5 * cex.lab)
    mtext(xlab, side = 3, line = -2, cex = cex.lab)
    
    par(mar = mar) 
    if (dr == "fert") {
      agev <- x$fert$functs$age
      ydr <- x$fert$functs$fert
    } else {
      agev <- x$surv$functs$age
      ydr <- x$surv$functs[[dr]]
    }
    nage <- length(agev)
    if ("xlim" %in% names(args)) {
      xlim <- args$xlim
    } else {
      xlim <- range(agev) * c(0, 1.1)
    }
    if ("ylim" %in% names(args)) {
      ylim <- args$ylim
    } else {
      if (dr == "surv") {
        ylim <- c(0, 1)
      } else {
        ylim <- c(0, max(ydr))
      }
    }
    plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
    
    lines(agev, ydr, type = type, col = col, lwd = lwd)
    Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
    Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
    mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  }
}

# ----------------------- #
# LIFE HISTORY VARIABLES: 
# ----------------------- #
# Main life history function:
CalcLifeHist <- function(theta = NULL, beta = NULL, dx = NULL, 
                         model = "GO", shape = "simple",
                         modelFert = "quadratic", ageMatur = 0, 
                         maxAge = NULL, lambdaMethod = "matrix") {
  
  # ----------------------------------------------- #
  # CHECK THAT THE MINIMUM INFORMATION IS PROVIDED:
  # ----------------------------------------------- #
  if (all(is.null(theta), is.null(beta))) {
    stop("No basic information provided. 
    Provide vectors of survival and fertility parameters 'theta' and 'beta'.")
  }
  
  # ------------------------------------------------ #
  # A) GENERAL SETUP AND VERIFICATION OF PARAMETERS:
  # ------------------------------------------------ #
  # Verify fertility model:
  .VerifyFertMod(modelFert = modelFert)
  
  # Extract beta attributes:
  betaAttr <- .SetBeta(beta = beta, modelFert = modelFert)
  beta <- betaAttr$beta
  
  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)
  
  # Extract theta attributes:
  thetaAttr <- .SetTheta(theta = theta, model = model, shape = shape)
  theta <- thetaAttr$theta
  
  # --------------------- #
  # B) PREPARE FUNCTIONS:
  # --------------------- #
  # a) Mortality function:
  .CalcMort <- function(theta, ...) UseMethod(".CalcMort")
  .CalcMort.matrix <- .DefineMortMatrix(model = model, shape = shape)
  .CalcMort.numeric <- .DefineMortNumeric(model = model, shape = shape)
  
  # b) Cummulative hazard:
  .CalcCumHaz <- function(theta, ...) UseMethod(".CalcCumHaz")
  .CalcCumHaz.matrix <- .DefineCumHazMatrix(model = model, shape = shape)
  .CalcCumHaz.numeric <- .DefineCumHazNumeric(model = model, shape = shape)
  
  # c) Survival:
  .CalcSurv <- function(theta, x) {
    exp(-.CalcCumHaz(theta = theta, x = x))
  }
  
  # d) Fertility function:
  .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
  .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)
  .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
  
  # ------------------------ #
  # PREPARE AGE INFORMATION:
  # ------------------------ #
  # age interval length:
  if (is.null(dx)) {
    dx <- 1
  }
  
  # Maximum age:
  if (is.null(maxAge)) {
    maxAge <- .FindMaxAge(theta, bound = 1e-10, survFunc = .CalcSurv)
  }
  
  #  age vector:
  x <- seq(0, maxAge, dx)
  
  # Number of age intervals:
  nx <- length(x)
  
  # Calcualte survival:
  Sx <- CalcSurv(theta = theta, x = x, model = model, shape = shape)
  
  # Age-specific survival:
  px <- CalcSurv(theta = theta, x = x + 1, model = model, shape = shape) /
    CalcSurv(theta = theta, x = x, model = model, shape = shape)
  
  # Age-specific fertility:
  mx <- c(rep(0, length(which(x < ageMatur))), 
          CalcFert(beta = beta, x = x[x >= ageMatur] - ageMatur, 
                   modelFert = modelFert))
  
  # 'Continous' age interval length:
  dxc <- 0.001
  
  # 'Continous' age vector:
  xc <- seq(0, maxAge, dxc)
  
  # 'Continuous' survival:
  Sxc <- CalcSurv(theta = theta, x = xc, model = model, shape = shape)
  
  # 'Continous' fertility:
  mxc <- c(rep(0, length(which(xc < ageMatur))), 
           CalcFert(beta = beta, x = xc[xc >= ageMatur] - ageMatur, 
                    modelFert = modelFert))
  
  # ------------------------------------------------ #
  # FIND LAMBDA, AGE STRUCTURE, REPRODUCTIVE VALUES:
  # ------------------------------------------------ #
  if (lambdaMethod == "Lotka") {
    # ------------------------------------- #
    # Using Lotka's (1913) renewal equation:
    # ------------------------------------- #
    Findr <- function(r) {
      lhs <- sum(exp(-r * xc) * Sxc * mxc * dxc)
      return((lhs - 1)^2)
    }
    
    # Find r:
    rout <- optimize(f = Findr, interval = c(0, 5), tol = 10e-20)
    
    # Extract r:
    r <- rout$minimum
    
    # Population growth rate:
    lambda <- exp(r)
    
    # Calculate stable age structure:
    w <- (Sx * exp(-r * x)) / sum(Sx * exp(-r * x))
    
    # Reproductive value:
    v <- sapply(x, function(xx) {
      idsx <- which(xc >= xx)
      sum(exp(-r * (xc[idsx] - xx)) * Sxc[idsx] * mxc[idsx] * dxc) / 
        Sxc[idsx[1]]
    })
    
  } else if (lambdaMethod == "matrix") {
    # ---------------------- #
    # Using matrix algebra:
    # ---------------------- #
    # Transition matrix:
    A <- .FillLeslieMatr(p = px, f = mx, n = nx)
    
    # Eigen analysis:
    eA <- eigen(A)
    eAt <- eigen(t(A))
    
    # Population growth rate:
    lambda <- Re(eA$value[1])
    w <- abs(Re(eA$vector[, 1])) / sum(abs(Re(eA$vector[, 1])))
    v <- abs(Re(eAt$vectors[, 1]))
    v <- v / v[1]
    r <- log(lambda)
  } else {
    stop("Wrong 'lambdaMethod', options are 'Lotka' or 'matrix'.")
  }
  
  # Net reproductive rate:
  R0 <- sum(Sxc * mxc * dxc)
  
  # Cohort generation time:
  Tc <- sum(xc * Sxc * mxc * dxc) / R0
  
  # Demographic dispersion:
  sigmad <- sum((xc - Tc)^2 * Sxc * mxc * dxc) / R0
  
  # Generation time of stable population:
  Ts <- sum(exp(-r * xc) * xc * Sxc * mxc * dxc)
  
  # Damping ratio:
  # tau <- lambda / abs(Re(eA$values[2]))
  
  # Population entropy:
  phi <- Sxc * mxc * lambda^(-xc) 
  idn0 <- which(phi > 0)
  phi <- phi[idn0]
  H <- - 1 / Ts * sum(phi * log(phi) * dxc)
  
  # Sensitivities:
  # spx <- v * w / sum(v * w)
  # smx <- v[1] * w / sum(v * w)
  spx <- c(v[-1] * w[-nx] / sum(v * w), NA)
  smx <- c(v[1] * w[-nx] / sum(v * w), NA)
  smx[which(x < ageMatur)] <- 0
  
  # Elasticities:
  epx <- px / lambda * spx
  emx <- mx / lambda * smx
  
  # Prepare outputs:
  lhlabs <- c(r = "Intrinsic population growth rate", 
              lambda = "Stable population growth", 
              Ts = "Generation time of stable population", 
              Tc = "Cohort generation time", 
              R0 = "Net reproductive rate", 
              sigmad = "Demographic dispersion", 
              H = "Population entropy")
  lhnames <- names(lhlabs)
  
  lifeHist <- data.frame(Description = lhlabs, Variable = lhnames, 
                         Value = c(r, lambda, Ts, Tc, R0, sigmad, H))
  rownames(lifeHist) <- NULL
  
  # Start output list:
  outList <- list(demoTab = data.frame(Age = x, Sx = Sx, px = px, qx = 1 - px,
                                       mx = mx, w = w, v = v, sp = spx, 
                                       sm = smx, ep = epx, em = emx),
                  lifeHist = lifeHist,
                  settings = list(theta = theta, beta = beta, model = model,
                                  shape = shape, modelFert = modelFert,
                                  ageMatur = ageMatur))
  # Assign class:
  class(outList) <- c("PDlifeHist")
  
  # return output:
  return(outList)
}

# Print life history:
print.PDlifeHist <- function(x, ...) {
  demolabs <- c(px = "Age-specific survival", 
                mx = "Age-specific fertility", 
                w = "Stable age structure",
                v = "Reproductive value",
                sp = "Sensitivity to survival",
                sm = "Sensitivity to fertility",
                ep = "Elasticity to survival",
                em = "Elasticity to fertility")
  
  
  cat("Life history variables:\n")
  cat("-----------------------\n")
  print(x$lifeHist, ...)
  
}

# Summary life history:
summary.PDlifeHist <- function(object, ...) {
  
  cat("Settings:\n")
  cat("=========\n")
  cat("Survival:\n")
  cat("---------\n")
  cat(sprintf("model: %s\n", object$settings$model))
  cat(sprintf("shape: %s\n", object$settings$shape))
  cat("theta parameters:\n")
  print(object$settings$theta)
  cat("\nFertility:\n")
  cat("----------\n")
  cat(sprintf("modelFert: %s\n", object$settings$modelFert))
  cat("beta parameters:\n")
  print(object$settings$beta)
  
  
  cat("\nLife history variables:\n")
  cat("=======================\n")
  print(object$lifeHist, ...)
  
}

# plot life history:
plot.PDlifeHist <- function(x, type = "rates", ...) {
  # User par settings:
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Additional arguments:
  args <- list(...)
  nameArgs <- names(args)
  
  # demo labels:
  demolabs <- c(px = expression(paste("Age-specific survival, ", 
                                      italic(p[x]))), 
                mx = expression(paste("Age-specific fertility, ", 
                                      italic(m[x]))), 
                w = expression(paste("Stable age structure, ", omega)),
                v = expression(paste("Reproductive value, ", nu)),
                sp = expression(paste("Sensitivity to survival, ", 
                                      italic(s[p]))),
                sm = expression(paste("Sensitivity to fertility, ", 
                                      italic(s[m]))),
                ep = expression(paste("Elasticity to survival, ", 
                                      italic(e[p]))),
                em = expression(paste("Elasticity to fertility, ", 
                                      italic(e[m]))))
  
  # xlim:
  if ("xlim" %in% nameArgs) {
    xlim <- args$xlim
  } else {
    xlim <- range(x$demoTab$Age)
  }
  
  # Color:
  if ("col" %in% nameArgs) {
    col <- args$col 
  } else {
    col <- "dark red"
  }
  
  # lwd:
  if ("lwd" %in% nameArgs) {
    lwd <- args$lwd
  } else {
    lwd <- 2
  }
  if (type == "rates") {
    plList <- c("px", "mx", "v", "w")
  } else if (type == "sensitivity") {
    plList <- c("sp", "sm", "ep", "em")
  } else if (type == "all") {
    plList <- c("px", "mx", "w", "v", "sp", "sm", "ep", "em")
  } else {
    stop("Wrong 'type' argument, options are 'rates', 'sensitivity', or 'all'.")
  }
  npl <- length(plList)
  
  # Demographic rates:
  # ------------------ #
  # ylim:
  if ("ylim" %in% nameArgs) {
    ylim <- args$ylim
  } else {
    ylim <- sapply(plList, function(ix) {
      if (ix %in% c("sp", "sm")) {
        yl <- c(0, max(c(x$demoTab$sp, x$demoTab$sm), na.rm = TRUE))
      } else if (ix %in% c("ep", "em")) {
        yl <- c(0, max(c(x$demoTab$ep, x$demoTab$em), na.rm = TRUE))
      } else {
        yl <- c(0, max(x$demoTab[[ix]]))
      }
      return(yl)
    })
  }
  
  par(mfrow = c(npl / 2, 2), mar = c(4.1, 4.1, 1, 1))
  for (demi in plList) {
    if (inherits(ylim, "matrix")) {
      yylim <- ylim[, demi]
    } else {
      yylim <- ylim
    }
    
    plot(x$demoTab$Age, x$demoTab[[demi]], type = 'l', col = col, lwd = lwd, 
         xlim = xlim, ylim = yylim, xlab = "Age", ylab = demolabs[demi])
  }
  
}

# --------------------------------------------------- #
# SUMMARY STATISTICS AND OTHER DEMOGRAPHIC VARIABLES:
# --------------------------------------------------- #
# Function to calculate remaining life expectancy:
CalcRemainLifeExp <- function(theta, x = NULL, dx = NULL, xmax = NULL,
                              atAllAges = FALSE, model = "GO", 
                              shape = "simple", checkTheta = TRUE) {
  # ------------------------------------------------
  # A) GENERAL SETUP AND VERIFICATION OF PARAMETERS:
  # ------------------------------------------------
  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)

  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }

  # ---------------------
  # B) PREPARE FUNCTIONS:
  # ---------------------
  # b) Cummulative hazard:
  .CalcCumHaz <- function(theta, ...) UseMethod(".CalcCumHaz")
  .CalcCumHaz.matrix <- .DefineCumHazMatrix(model = model, shape = shape)
  .CalcCumHaz.numeric <- .DefineCumHazNumeric(model = model, shape = shape)

  # c) Survival:
  .CalcSurv <- function(theta, x) {
    exp(-.CalcCumHaz(theta = theta, x = x))
  }

  # Calculate remaining life expectancy:
  if (is.null(xmax)) {
    if (is.matrix(theta)) {
      xmax <- max(apply(theta, 1, function(th) {
        .FindMaxAge(th, bound = 0.00001, survFunc = .CalcSurv)
      }))
    } else {
      xmax <- .FindMaxAge(theta, bound = 0.00001, survFunc = .CalcSurv)
    }
  }
  dx <- ifelse(is.null(dx), xmax / 10000, dx)
  if (is.null(x)) {
    x <- 0
  }
  xv <- seq(min(x), xmax, dx)
  nx <- length(xv)
  idx <- sapply(x, function(xx) {
    which(abs(xx - xv) == min(abs(xx - xv)))[1]
  })

  if (is.matrix(theta)) {
    Sx <- t(apply(theta, 1, function(th) {
      .CalcSurv(th, xv)
    }))
    Ex <- apply(Sx, 1, function(sx) {
      rev(cumsum(rev(sx))) * dx / sx
    })
    colnames(Ex) <- sprintf("RemLExpTheta%s", 1:nrow(theta))
    remLexp <- cbind(Age = xv, Ex)
  } else {
    Sx <- .CalcSurv(theta, xv)
    Ex <- rev(cumsum(rev(Sx * dx))) / Sx
    remLexp <- cbind(Age = xv, RemLExp = Ex)
  }
  if (!atAllAges) {
    remLexp <- remLexp[idx, ]
  }
  return(remLexp)
}

# Actuarial ageing rates:
CalcAgeingRateMort <- function(theta, x, model = "GO", shape = "simple",
                               checkTheta = TRUE) {
  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)
  
  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }
  
  # Ageing rate function:
  arFun <- .DefineARmort(model = model, shape = shape)
  
  # Calculate ageing rate:
  ar <- cbind(Age = x, AR = arFun(theta = theta, x = x))
  rownames(ar) <- NULL
  
  return(ar)
}

# Reproductive ageing rate:
CalcAgeingRateFert <- function(beta, x, modelFert = "quadratic", ageMatur = 0,
                               checkBeta = TRUE) {
  # Verify fertility model:
  .VerifyFertMod(modelFert = modelFert)
  
  # Extract beta attributes:
  if (checkBeta) {
    betaAttr <- .SetBeta(beta = beta, modelFert = modelFert)
    beta <- betaAttr$beta
  }
  
  # Verify that the age x is larger than the age at maturity:
  if (all(x < ageMatur)) {
    stop("Ages 'x' occur before the age at maturity 'ageMatur'.")
  } else {
    idlow <- which(x > ageMatur)
  }
  # a) Fertility method:
  .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
  .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
  .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)
  
  # Age increase:
  dx <- 0.000001
  
  # Fertility at x and x+dx
  fertx <- .CalcFert(beta = beta, x = x[idlow] - ageMatur)
  fertxdx <- .CalcFert(beta = beta, x = x[idlow] - ageMatur + dx)
  
  # Log-fertility:
  lnfx <- log(fertx)
  lnfxdx <- log(fertxdx)
  
  # Ageing rate:
  Ar <- rep(NA, length(x))
  Ar[idlow] <- (lnfxdx - lnfx) / dx
  fAr <- cbind(x = x, AR = Ar)
  return(fAr)
}

# Calculate age at maximum Fertility:
CalcAgeMaxFert <- function(beta, modelFert = "quadratic", ageMatur = 0, 
                           maxAge = 100) {
  
  # Verify fertility model:
  .VerifyFertMod(modelFert = modelFert)
  
  # Extract beta attributes:
  betaAttr <- .SetBeta(beta = beta, modelFert = modelFert)
  beta <- betaAttr$beta
  
  # a) Fertility method:
  .CalcFert <- function(beta, ...) UseMethod(".CalcFert")
  .CalcFert.matrix <- .DefineFertilityMatrix(modelFert = modelFert)
  .CalcFert.numeric <- .DefineFertilityNumeric(modelFert = modelFert)
  
  # Find age at maximum Fertility:
  xv <- seq(0, maxAge - ageMatur, 0.0001)
  if (modelFert %in% c("quadratic", "PeristeraKostaki")) {
    xm <- beta["b2"]
    dd <- 0
    ii <- 0
  } else if (modelFert %in% c("ColcheroMuller", "Hadwiger")) {
    if (modelFert == "ColcheroMuller") {
      dfdx <- function(x, beta) {
        x^3 + x^2 * (2 - beta["b2"]) + x * (1 - 2 * beta["b2"]) +
          beta["b3"] / (2 * beta["b1"]) - beta["b2"]
      }
      dfdx2 <- function(x, beta) {
        3 * x^2 + 2 * x * (2 - beta["b2"]) + 1 - 2 * beta["b2"]
      }
    } else {
      dfdx <- function(x, beta) {
        3/2 * x^(-1) - beta["b1"]^2 * beta["b2"] * x^(-2) + 
          beta["b1"]^2 / beta["b2"]
      }
      dfdx2 <- function(x, beta) {
        -3 / 2 * x^(-2) + 2 * beta["b1"]^2 * beta["b2"] * x^(-3)
      }
    }

    id0 <- which(sign(dfdx(xv[-length(xv)], beta)) != 
                   sign(dfdx(xv[-1], beta)))
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
  } else if (modelFert == "gamma") {
    if (beta["b1"] > 1) {
      xm <- (beta["b1"] - 1) / beta["b2"]
    } else {
      xm <- 0
    }
    dd <- 0
    ii <- 0
  } else if (modelFert == "beta") {
    if (beta["b1"] > 1 & beta["b2"] > 1) {
      xm <- beta["b3"] + (beta["b1"] - 1) / (beta["b1"] + beta["b2"] - 2) * 
        (beta["b4"] - beta["b3"])
    } else if (beta["b1"] < 1 & beta["b2"] < 1) {
      xm <- beta[c("b3", "b4")]
    } else if ((beta["b1"] < 1 & beta["b2"] >= 1) | 
               (beta["b1"] == 1 & beta["b2"] > 1)) {
      xm <- beta["b3"]
    } else if ((beta["b1"] >= 1 & beta["b2"] < 1) | 
               (beta["b1"] > 1 & beta["b2"] == 1)) {
      xm <- beta["b4"]
    } else {
      xm <- NA
    }
    dd <- 0
    ii <- 0
  }
  maxFert <- .CalcFert(beta, xm)
  maxFertv <- c(xm + ageMatur, maxFert, dd, ii, ageMatur)
  names(maxFertv) <- c("Age", "maxFert", "error", "iterations", "ageMatur")
  return(maxFertv)
}


# ---------------------------------------- #
# DISCRETE AGE OR STAGE DEMOGRAPHIC RATES:
# ---------------------------------------- #
# Demographic rates in discrete age intervals:
CalcDiscrDemo <- function(demo, dx = 1) {
  if (!(inherits(demo, "demoSurv") | inherits(demo, "demoBoth"))) {
    stop("Object 'demo' should be of class 'demoSurv'.\nCreate object demo with function CalcDemo().", call. = FALSE)
  }
  if (length(demo$surv$functs$age) == 1) {
    stop("Age-specific probabilities cannot be calculated from a single age.\nIncrease the age vector on CalcDemo().", call. = FALSE)
    
  }
  # ---------------------------- #
  # ---- Survival probabilities:
  # ---------------------------- #
  x <- demo$surv$functs$age
  mux <- demo$surv$functs$mort
  cumh <- demo$surv$functs$cumhaz
  Sx <- demo$surv$functs$surv
  xdis <- x[which(x %in% 0:max(x))]
  idAges <- which(x %in% xdis)
  nAges <- length(idAges)
  px <- Sx[idAges[-1]] / Sx[idAges[-nAges]]
  qx <- 1 - px
  lx <- Sx[idAges[-nAges]]
  demoDiscr <- cbind(age = xdis[-nAges], lx, px, qx)
  
  # ---------- #
  # Fertility:
  # ---------- #
  if (inherits(demo, "demoBoth")) {
    idMidAges <- which(c(x %in% (xdis + dx / 2)))
    if (length(idMidAges) == 0) {
      idMidAges <- idAges * 0
      for (ai in 1:nAges) {
        xdmid <- xdis + dx / 2
        idd <- which(abs(x - xdmid) == min(x - xdmid))
        idMidAges[ai] <- idd
      }
    }
    bx <- demo$fert$functs$fert[idMidAges[-nAges]] * dx
    demoDiscr <- cbind(demoDiscr, bx = bx)
  }
  class(demoDiscr) <- c("discrDemo", "matrix")
  return(demoDiscr)
}

# --------------------- #
# LIFE TABLE FUNCTIONS: 
# --------------------- #
# Calculate life table:
CalcLifeTable <- function(ageLast, ageFirst = NULL, departType, dx = 1, 
                          calcCIs = FALSE, nboot = 1000, alpha = 0.05) {
  # Error handling:
  if (missing(ageLast)) {
    stop("Argument 'ageLast' missing. 
         Provide vector of ages at last observation.")
  }
  
  # Number of records:
  n <- length(ageLast)
  
  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, n)
  }
  
  if (missing(departType)) {
    stop("Argument 'departType' missing. 
         Provide character vector of departure types (i.e., 'D' and 'C')")
  }
  
  if (n != length(departType)) {
    stop("Lengths of 'ageLast' and 'departType' differ.")
  } else if (n != length(ageFirst)) {
    stop("Lengths of 'ageLast' and 'ageFirst' differ.")
  }
  
  # Produce life table:
  ltMean <- .CalcLT(ageLast = ageLast, ageFirst = ageFirst, 
                    departType = departType, dx = dx)
  
  # Check for negative survivals:
  if (any(ltMean[, "lx"] < 0)) {
    warning("Negative survivals produced due to excess truncation in some age intervals.
            Consider better using function 'CalcProductLimitEst'")
  }
  
  if (calcCIs) {
    # Number of ages:
    nage <- nrow(ltMean)
    
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
      ltb <- .CalcLT(ageLast = ageLastBoot, ageFirst = ageFirstBoot, 
                     departType = departTypeBoot, dx = dx)
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
    ltcis$Settings <- c(calc = 1, nboot = nboot, alpha = alpha)
    
  } else {
    ltcis <- list(lx = NA, qx = NA, px = NA, ex = NA, 
                  Settings = c(calc = 0, nboot = NA, alpha = NA))
  }
  lt <- list(lt = ltMean, CIs = ltcis)
  class(lt) <- c("paramDemoLT")
  return(lt)
}

# Plot lifetable:
plot.paramDemoLT <- function(x, demorate = "lx", inclCIs = FALSE, ...) {
  # User par settings:
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Additional arguments:
  args <- list(...)
  
  # Labels for demographic rates:
  demRates <- c("lx", "qx", "px", "ex")
  
  # Check if demorate is properly provided:
  if (!demorate %in% c(demRates, "all")) {
    stop("Argument 'demorate' incorrect, values should be 'lx',\n'qx', 'px', 'ex', or 'all'.", call. = FALSE)
  }
  
  # Plot names for demographic rates:
  demRatesNames <- c(lx = expression(paste("Survival, ", italic(l[x]))),
                     qx = expression(paste("Mortality probability, ", 
                                           italic(q[x]))), 
                     px = expression(paste("Survival probability, ", 
                                           italic(p[x]))),
                     ex = expression(paste("Remaining life expectancy, ", 
                                           italic(e[x]))))
  
  # mar values:
  if ("mar" %in% names(args)) {
    mar <- args$mar
  } else {
    mar <- c(2, 4, 1, 1)
  }
  
  # Axis labels:
  if ("xlab" %in% names(args)) {
    xlab <- args$xlab
  } else {
    xlab <- expression(paste("Age, ", italic(x)))
  }
  if ("ylab" %in% names(args)) {
    ylab <- args$ylab
    if (demorate == "all" & length(ylab) != 4) {
      stop("Length of ylab needs to be four for demorate = 'all'.\nNote variables plotted will be 'lx', 'qx', 'px', and 'ex'.", 
           call. = FALSE)
    }
  } else {
    if (demorate == "all") {
      ylab <- demRatesNames
    } else {
      ylab <- demRatesNames[demorate]
    }
  }
  
  # Type of line:
  if ("type" %in% names(args)) {
    type <- args$type
  } else {
    type <- "l"
  }
  
  # Color
  if ("col" %in% names(args)) {
    col <- args$col
  } else {
    col <- "#800026"
  }
  
  # Line width:
  if ("lwd" %in% names(args)) {
    lwd <- args$lwd
  } else {
    lwd <- 1
  }
  
  # Cex.lab:
  if ("cex.lab" %in% names(args)) {
    cex.lab <- args$cex.lab
  } else {
    cex.lab <- 1
  }
  
  # Cex.axis:
  if ("cex.axis" %in% names(args)) {
    cex.axis <- args$cex.axis
  } else {
    cex.axis <- 1
  }
  
  if (demorate == "all") {
    # -------------- #
    # MULTIPLE PLOTS:
    # -------------- #
    # Layout matrix:
    laymat <- cbind(c(2, 4, 1), c(3, 5, 1))
    
    # Widths and heights:
    widths <- c(1, 1)
    heights <- c(1, 1, 0.15)
    
    # produce plot:
    layout(mat = laymat, widths = widths, heights = heights)
    
    # x-axis label:
    par(mar = mar * c(0, 1, 0, 1))
    plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    text(0.5, 0.5, xlab, cex = 1.5 * cex.lab)
    
    par(mar = mar) 
    for (dr in demRates) {
      agev <- x$lt[, "Ages"]
      nage <- length(agev)
      if ("xlim" %in% names(args)) {
        xlim <- args$xlim
      } else {
        xlim <- range(agev) * c(0, 1.1)
      }
      ydr <- x$lt[, dr]
      if ("ylim" %in% names(args)) {
        ylim <- args$ylim
      } else {
        if (dr != "ex") {
          ylim <- c(0, 1)
        } else {
          ylim <- c(0, max(ydr))
        }
      }
      if (dr == "lx") {
        type <- "s"
      } else {
        type <- "l"
      }
      plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
      if (inclCIs & x$CIs$Settings["calc"] == 1) {
        ydrCI <- x$CIs[[dr]]
        colci <- adjustcolor(col, alpha.f = 0.25)
        if (dr == "lx") {
          for (ii in 1:nage) {
            xp <- agev[ii] + c(0, 1)
            yp <- ydrCI[ii, c("Lower", "Upper")]
            polygon(xp[c(1, 2, 2, 1)], yp[c(1, 1, 2, 2)], col = colci, 
                    border = NA)
          }
        } else {
          xp <- agev
          yp <- ydrCI[, c("Lower", "Upper")]
          polygon(c(xp, rev(xp)), c(yp[, 1], rev(yp[, 2])), col = colci, 
                  border = NA)
        }
        
      }
      
      lines(agev, ydr, type = type, col = col, lwd = lwd)
      Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
      Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
      mtext(ylab[dr], side = 2, line = mar[2] / 2, cex = cex.lab)
    }    
  } else {
    # ------------ #
    # SINGLE PLOT:
    # ------------ #
    dr <- demorate
    # Layout matrix:
    laymat <- matrix(c(2, 1), nrow = 2, ncol = 1)
    
    # Widths and heights:
    widths <- c(1)
    heights <- c(1, 0.15)
    
    # produce plot:
    layout(mat = laymat, widths = widths, heights = heights)
    
    # x-axis label:
    par(mar = mar * c(0, 1, 0, 1))
    plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
    # text(0.5, 0.5, xlab, cex = 1.5 * cex.lab)
    mtext(xlab, side = 3, line = -2, cex = cex.lab)
    
    par(mar = mar) 
    agev <- x$lt[, "Ages"]
    nage <- length(agev)
    if ("xlim" %in% names(args)) {
      xlim <- args$xlim
    } else {
      xlim <- range(agev) * c(0, 1.1)
    }
    ydr <- x$lt[, dr]
    if ("ylim" %in% names(args)) {
      ylim <- args$ylim
    } else {
      if (dr != "ex") {
        ylim <- c(0, 1)
      } else {
        ylim <- c(0, max(ydr))
      }
    }
    if (dr == "lx") {
      type <- "s"
    } else {
      type <- "l"
    }
    plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
    if (inclCIs & x$CIs$Settings["calc"] == 1) {
      ydrCI <- x$CIs[[dr]]
      colci <- adjustcolor(col, alpha.f = 0.25)
      if (dr == "lx") {
        for (ii in 1:nage) {
          xp <- agev[ii] + c(0, 1)
          yp <- ydrCI[ii, c("Lower", "Upper")]
          polygon(xp[c(1, 2, 2, 1)], yp[c(1, 1, 2, 2)], col = colci, 
                  border = NA)
        }
      } else {
        xp <- agev
        yp <- ydrCI[, c("Lower", "Upper")]
        polygon(c(xp, rev(xp)), c(yp[, 1], rev(yp[, 2])), col = colci, 
                border = NA)
      }
      
    }
    lines(agev, ydr, type = type, col = col, lwd = lwd)
    Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
    Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
    mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  }
}

# ------------------------ #
# PRODUCT-LIMIT ESTIMATOR: 
# ------------------------ #
# Simple Product-limit estimator function:
CalcProductLimitEst <- function(ageLast, ageFirst = NULL, departType,
                                calcCIs = FALSE, nboot = 1000, alpha = 0.05) {
  # Error handling:
  if (missing(ageLast)) {
    stop("Argument 'ageLast' missing. 
         Provide vector of ages at last observation.")
  }
  
  # Number of records:
  n <- length(ageLast)
  
  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, n)
  } 
  
  # Eliminate potential rounding errors:
  ageFirst <- round(ageFirst, 8)
  ageLast <- round(ageLast, 8)
  
  if (missing(departType)) {
    stop("Argument 'departType' missing. 
         Provide character vector of departure types (i.e., 'D' and 'C')")
  }
  
  if (n != length(departType)) {
    stop("Lengths of 'ageLast' and 'departType' differ.")
  } else if (n != length(ageFirst)) {
    stop("Lengths of 'ageLast' and 'ageFirst' differ.")
  }
  
  # Find records with same first and last age:
  idsame <- which(ageLast == ageFirst)
  
  # Increase last age by one day:
  ageLast[idsame] <- ageLast[idsame] + 1/365.25
  
  # Product limit estimator:
  pleTab <- .CalcPLE(ageLast = ageLast, ageFirst = ageFirst,
                     departType = departType)
  
  # Confidence intervals:
  if (calcCIs) {
    # Number of individual ages:
    nple <- length(pleTab$Ages)
    pleAges <- pleTab$Ages
    
    # ==================== #
    # ==== BOOTSTRAP: ====
    # ==================== #
    bootTab <- matrix(NA, nple, nboot)
    for (iboot in 1:nboot) {
      idboot <- sample(1:n, size = n, replace = TRUE)
      ageLastb <- ageLast[idboot]
      ageFirstb <- ageFirst[idboot]
      departTypeb <- departType[idboot]
      
      pleBooti <- .CalcPLE(ageLast = ageLastb, ageFirst = ageFirstb,
                           departType = departTypeb)
      
      idages <- which(pleAges %in% pleBooti$Ages)
      bootTab[idages, iboot] <- pleBooti$ple
      idna <- which(is.na(bootTab[, iboot]))
      nna <- length(idna)
      while(nna > 0) {
        bootTab[idna, iboot] <- bootTab[idna - 1, iboot]
        idna <- which(is.na(bootTab[, iboot]))
        nna <- length(idna)
        
      }
    }
    
    # Calculate CIs:  # Calculate CIs:
    pleci <- t(apply(bootTab, 1, quantile, c(alpha / 2, 1 - alpha / 2),
                     na.rm = TRUE))
    colnames(pleci) <- c("Lower", "Upper")
    
    # Final data frame:
    pleTab <- data.frame(pleTab, pleci)
  }
  
  class(pleTab) <- c("paramDemoPLE", "data.frame")
  return(pleTab)
}


# Plot Product-limit estimator:
plot.paramDemoPLE <- function(x, inclCIs = FALSE, ...) {
  # User par settings:
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Additional arguments:
  args <- list(...)
  
  # Vector of ages:
  agev <- x$Ages
  nage <- length(agev)
  
  # mar values:
  if ("mar" %in% names(args)) {
    mar <- args$mar
  } else {
    mar <- c(2, 4, 1, 1)
  }
  
  # Colors:
  if ("col" %in% names(args)) {
    col <- args$col
  } else {
    col <- "#800026"
  }
  
  # Axis labels:
  if ("xlab" %in% names(args)) {
    xlab <- args$xlab
  } else {
    xlab <- expression(paste("Age, ", italic(x)))
  }
  if ("ylab" %in% names(args)) {
    ylab <- args$ylab
  } else {
    ylab <- "Product-limit estimator"
  }
  
  # Type of line:
  if ("type" %in% names(args)) {
    type <- args$type
  } else {
    type <- "l"
  }
  
  # Color
  if ("col" %in% names(args)) {
    col <- args$col
  } else {
    col <- "#800026"
  }
  
  # Line width:
  if ("lwd" %in% names(args)) {
    lwd <- args$lwd
  } else {
    lwd <- 1
  }
  
  # Cex.lab:
  if ("cex.lab" %in% names(args)) {
    cex.lab <- args$cex.lab
  } else {
    cex.lab <- 1
  }
  
  # Cex.axis:
  if ("cex.axis" %in% names(args)) {
    cex.axis <- args$cex.axis
  } else {
    cex.axis <- 1
  }
  
  if ("xlim" %in% names(args)) {
    xlim <- args$xlim
  } else {
    xlim <- range(agev) * c(0, 1.1)
  }
  if ("ylim" %in% names(args)) {
    ylim <- args$ylim
  } else {
    ylim <- c(0, 1)
  }
  
  # Logical in case CIs were calculated:
  CIsCalc <- ifelse ("Lower" %in% colnames(x), TRUE, FALSE)
  
  # Produce plots:
  par(mar = c(4, 4, 1, 1), mfrow = c(1, 1))
  plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
  if (inclCIs & CIsCalc) {
    for (ici in c("Lower", "Upper")) {
      lines(agev, x[[ici]], type = type, col = col, lwd = lwd, lty = 2)
    }
  }
  lines(agev, x$ple, type = type, col = col, lwd = lwd)
  Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
  Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
  mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  mtext(xlab, side = 1, line = 2, cex = cex.lab)
}

# --------------------------------------- #
# SAMPLING AND OTHER ANCILLARY FUNCTIONS:
# --------------------------------------- #
# Sampling random ages at death:
SampleRandAge <- function(n, theta, dx = 0.001, model = "GO", shape = "simple",
                          minSx = 0.0001) {
  
  # Extract demographic functions:
  demo <- CalcDemo(theta, dx = dx, model = model, shape = shape, minSx = minSx,
                   summarStats = FALSE, type = "survival")
  
  # Extract CDF of ages at death:
  Fx <- 1 - demo$surv$functs$surv
  
  # Draw random uniform values:
  u <- runif(n)
  
  # Extract ages for Fx = u:
  uages <- demo$surv$functs$age[findInterval(u, Fx, rightmost.closed = TRUE)]
  
  # return random ages:
  return(uages)
}

# ================================== END ===================================== #