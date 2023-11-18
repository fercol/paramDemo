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
    BETUPP <- all(sapply(1:nBe, function(bi) {
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

# ------------- #
# AGEING RATES:
# ------------- #
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

# --------------------------------------------- #
# EXTRACT SILER PARAMETERS FROM SURVIVAL RATES: 
# --------------------------------------------- #
.FindSilerPars <- function(theta, palpha) {
  theta[c(2, 3, 5)] <- abs(theta[c(2, 3, 5)])
  
  # Extract mortality, survival, etc:
  demotest <- CalcDemo(theta = theta, shape = "bathtub", minSx = 0.00001, 
                       type = "survival", summarStats = FALSE)
  
  # Vector of demographic variables:
  paltest <- c(CalcAveDemo(demotest), omega = max(demotest$age))
  
  # max alpha:
  alphaMax <- max(paltest[3], palpha[3])
  
  # Scale alpha's:
  paltest[3] <- paltest[3] / alphaMax
  palpha[3] <- palpha[3] / alphaMax
  
  # max maxAge:
  omegaMax <- max(paltest[4], palpha[4])
  
  # Scale max age:
  paltest[4] <- paltest[4] / omegaMax
  palpha[4] <- palpha[4] / omegaMax
  
  # Difference between
  diffpa <- palpha - paltest
  
  return(sum(diffpa^2))
}

# ---------------- #
# RECURSIVE OPTIM:
# ---------------- #
.myoptim <- function(par, fn, ...) {
  opt1 <- optim(par = par, fn = fn, ...)
  ntry <- 0
  while(opt1$convergence != 0 & ntry < 50) {
    ntry <- ntry + 1
    opt1 <- optim(par = opt1$par, fn = fn, ...)
  }
  return(opt1)
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
    } else {
      ageMaxFert <- NA
      calc <- FALSE
    }
    fertList <- list(functs = data.frame(age = x[which(x >= ageMatur)], 
                                         fert = fert),
                     summStats = list(calculated = calc, 
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
  op <- par(no.readonly = TRUE)
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
  op <- par(no.readonly = TRUE)
  
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
  par(op)
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
      dfdx <- function(x0, beta) {
        x0^3 + x0^2 * (2 - beta["b2"]) + x0 * (1 - 2 * beta["b2"]) +
          beta["b3"] / (2 * beta["b1"]) - beta["b2"]
      }
      dfdx2 <- function(x0, beta) {
        3 * x0^2 + 2 * x0 * (2 - beta["b2"]) + 1 - 2 * beta["b2"]
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
# Age-specific survival probability:
CalcSurvProbs <- function(demo, dx = 1) {
  if (!class(demo)[2] %in%  c("demoSurv", "demoBoth")) {
    stop("Object 'demo' should be of class 'demoSurv'.\nCreate object demo with function CalcDemo().", call. = FALSE)
  }
  if (length(demo$age) == 1) {
    stop("Age-specific probabilities cannot be calculated from a single age.\nIncrease the age vector on CalcDemo().", call. = FALSE)

  }
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
  demoProbs <- cbind(age = xdis[-nAges], px, qx, lx)
  class(demoProbs) <- "demoProbs"
  return(demoProbs)
}

# Extract average juvenile and adult vital rates:
CalcAveDemo <- function(demo) {
  
  # Find life table probabilities
  probs <- CalcSurvProbs(demo)
  
  # Find age at minimum mortality:
  alpha <- floor(demo$surv$functs$age[which(demo$surv$functs$mort == 
                                              min(demo$surv$functs$mort))])
  
  # Indices of adult and juvenile ages:
  idjuv <- which(probs[, "age"] < alpha)
  idad <- which(probs[, "age"] >= alpha)
  
  # Calculate average survival for juvs and ads:
  pj <- sum(probs[idjuv, "px"] * probs[idjuv, "lx"]) /
    sum(probs[idjuv, "lx"])
  pa <- sum(probs[idad, "px"] * probs[idad, "lx"]) / sum(probs[idad, "lx"])
  
  # Vector of demographic variables:
  palpha <- c(pj = pj, pa = pa, alpha = alpha)
  
  # return vector:
  return(palpha)
}

# --------------------- #
# LIFE TABLE FUNCTIONS: 
# --------------------- #
# Calculate life table:
CalcLifeTable <- function(ageLast, ageFirst = NULL, departType, dx = 1) {
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
    # Dx[ix] <- sum(lived[idDx])
    Dx[ix] <- nDx
    
    
    # C) PROPORTION LIVED BY THOSE THAT DIED IN INTERVAL:
    if (Dx[ix] > 1) {
      ylived <- xLast[idDx] - agev[ix]
      # ax[ix] <- sum(ylived) / Dx[ix]
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
  
  # Life-table:
  lt <- cbind(Ages = agev, Nx = Nx, Dx = Dx, lx = lx, px = px,
              qx = qx, Lx = Lx, Tx = Tx, ex = ex)
  class(lt) <- c("paramDemoLT")
  return(lt)
}

# Calculation of life table CIs:
CalcLifeTableCIs <- function(ageLast, ageFirst = NULL, departType, 
                             nboot = 2000, alpha = 0.05, dx = 1) {
  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, n)
  }
  
  # Unit age vector for that sex:
  agev <- seq(from = 0, to = ceiling(max(ageLast[which(departType == "D")])), 
              by = dx)
  nage <- length(agev)
  
  # Mean life table:
  ltMean <- CalcLifeTable(ageLast = ageLast, ageFirst = ageFirst, 
                          departType = departType, dx = dx)
  
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
  ltcis$Settings <- c(nboot = nboot, alpha = alpha)
  
  # class:
  class(ltcis) <- "paramDemoLTCIs"
  return(ltcis)
}

# Plot lifetable:
plot.paramDemoLT <- function(x, demorate = "lx", ...) {
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
  
  # Produce plots:
  op <- par(no.readonly = TRUE)
  
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
      agev <- x[, "Ages"]
      nage <- length(agev)
      if ("xlim" %in% names(args)) {
        xlim <- args$xlim
      } else {
        xlim <- range(agev) * c(0, 1.1)
      }
      ydr <- x[, dr]
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
    agev <- x[, "Ages"]
    nage <- length(agev)
    if ("xlim" %in% names(args)) {
      xlim <- args$xlim
    } else {
      xlim <- range(agev) * c(0, 1.1)
    }
    ydr <- x[, dr]
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
    
    lines(agev, ydr, type = type, col = col, lwd = lwd)
    Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
    Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
    mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  }
  par(op)
}

# Plotting life table CIs:
plot.paramDemoLTCIs <- function(x, demorate = "lx", ...) {
  
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
  
  # Line width:
  if ("lwd" %in% names(args)) {
    lwd <- args$lwd
  } else {
    lwd <- 1
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
  
  # Colors:
  if ("col" %in% names(args)) {
    if (length(args$col) < 2) {
      cols <- c(mean = "#800026", cis = "#FC4E2A")
    } else {
      cols <- args$col[1:2]
      names(cols) <- c("mean", "cis")
    }
  } else {
    cols <- c(mean = "#800026", cis = "#FC4E2A")
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
  
  # Cex.legend:
  if ("cex.legend" %in% names(args)) {
    cex.legend <- args$cex.legend
  } else {
    cex.legend <- 1
  }
  
  # Produce plots:
  op <- par(no.readonly = TRUE)
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
    legend('right', c("Life table", "CIs"), col = cols, pch = c(NA, 15), 
           lwd = c(lwd, NA), pt.cex = c(NA, 2), bty = "n", cex = cex.legend)
    
    par(mar = mar) 
    for (dr in demRates) {
      agev <- x[[dr]][, "Ages"]
      nage <- length(agev)
      if ("xlim" %in% names(args)) {
        xlim <- args$xlim
      } else {
        xlim <- range(agev) * c(0, 1.1)
      }
      ydr <- x[[dr]][, -1]
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
      if (dr == "lx") {
        for (ii in 1:nage) {
          xp <- agev[ii] + c(0, 1)
          yp <- ydr[ii, c("Lower", "Upper")]
          polygon(xp[c(1, 2, 2, 1)], yp[c(1, 1, 2, 2)], col = cols["cis"], 
                  border = NA)
        }
      } else {
        xp <- agev
        yp <- ydr[, c("Lower", "Upper")]
        polygon(c(xp, rev(xp)), c(yp[, 1], rev(yp[, 2])), col = cols["cis"], 
                border = NA)
      }
      lines(agev, ydr[, dr], type = type, col = cols["mean"], lwd = lwd)
      Axis(xlim + c(0, 10), side = 1, pos = ylim[1], cex.axis = cex.axis)
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
    legend('right', c("Life table", "CIs"), col = cols, pch = c(NA, 15), 
           lwd = c(lwd, NA), pt.cex = c(NA, 2), bty = "n", cex = cex.legend)
    
    par(mar = mar) 
    agev <- x[[dr]][, "Ages"]
    nage <- length(agev)
    if ("xlim" %in% names(args)) {
      xlim <- args$xlim
    } else {
      xlim <- range(agev) * c(0, 1.1)
    }
    ydr <- x[[dr]][, -1]
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
    if (dr == "lx") {
      for (ii in 1:nage) {
        xp <- agev[ii] + c(0, 1)
        yp <- ydr[ii, c("Lower", "Upper")]
        polygon(xp[c(1, 2, 2, 1)], yp[c(1, 1, 2, 2)], col = cols["cis"], 
                border = NA)
      }
    } else {
      xp <- agev
      yp <- ydr[, c("Lower", "Upper")]
      polygon(c(xp, rev(xp)), c(yp[, 1], rev(yp[, 2])), col = cols["cis"], 
              border = NA)
    }
    lines(agev, ydr[, dr], type = type, col = cols["mean"], lwd = lwd)
    Axis(xlim + c(0, 10), side = 1, pos = ylim[1], cex.axis = cex.axis)
    Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
    mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  }
  par(op)
  
}

# ------------------------ #
# PRODUCT-LIMIT ESTIMATOR: 
# ------------------------ #
# Simple Product-limit estimator function:
CalcProductLimitEst <- function(ageLast, ageFirst = NULL, departType) {
  # Number of individuals in dataset:
  n <- length(ageLast)
  
  # Find records with same first and last age:
  idsame <- which(ageLast == ageFirst)
  
  # Increase last age by one day:
  ageLast[idsame] <- ageLast[idsame] + 1/365.25
  
  # create identities and age list:
  if (is.null(ageFirst)) {
    allAges <- sort(ageLast)
    allAgesId <- departType
  } else {
    idAgeFirst <- which(ageFirst > min(ageLast))
    allAges <- c(ageFirst[idAgeFirst], ageLast)
    allAgesId <- c(rep("F", length(idAgeFirst)), departType)
  }
  
  allAgeTypes <- c("D", "C", "F")
  NallTypes <- length(allAgeTypes)
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
  recTab <- matrix(0, nages, NallTypes, dimnames = list(ageNames, allAgeTypes))
  for (at in ageTypes) {
    idtemp <- rep(0, nAllAges)
    idEqAt <- which(allAgesId == at)
    idtemp[idEqAt] <- 1
    ttemp <- table(allAges, idtemp)
    temp <- c(ttemp[, 2])
    recTab[rownames(ttemp), at] <- temp
  }
  
  # Cumulative tables:
  cumTab <- recTab * 0
  for (at in ageTypes) {
    cumTab[, at] <- rev(cumsum(rev(recTab[, at])))
  }
  
  Nx <- cumTab[, "D"] + cumTab[, "C"] - cumTab[, "F"]
  Dx <- recTab[, "D"]
  
  idDead <- which(recTab[, "D"] > 0)
  ple <- cumprod(1 - c(Dx / Nx)[idDead])
  
  # Add age 0:
  Ages <- unAllAges[idDead]
  if (Ages[1] > 0) {
    Ages <- c(0, Ages)
    ple <- c(1, ple)
  }
  
  # Fill up table:
  pleTab <- data.frame(Ages = Ages, ple = ple)
  class(pleTab) <- c("paramDemoPLE")
  return(pleTab)
}

# Product-limit estimator CIs:
CalcProductLimitEstCIs <- function(ageFirst, ageLast, departType, nboot = 1000,
                                 alpha = 0.05) {
  # =================== #
  # ==== NON-BOOT: ====
  # =================== #
  # Number of individuals in dataset:
  n <- length(ageLast)
  
  # Find records with same first and last age:
  idsame <- which(ageLast == ageFirst)
  
  # Increase last age by one day:
  ageLast[idsame] <- ageLast[idsame] + 1/365.25
  
  # create identities and age list:
  if (is.null(ageFirst)) {
    allAges <- sort(ageLast)
    allAgesId <- departType
  } else {
    idAgeFirst <- which(ageFirst > min(ageLast))
    allAges <- c(ageFirst[idAgeFirst], ageLast)
    allAgesId <- c(rep("F", length(idAgeFirst)), departType)
  }
  
  allAgeTypes <- c("D", "C", "F")
  NallTypes <- length(allAgeTypes)
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
  recTab <- matrix(0, nages, NallTypes, dimnames = list(ageNames, allAgeTypes))
  for (at in ageTypes) {
    idtemp <- rep(0, nAllAges)
    idEqAt <- which(allAgesId == at)
    idtemp[idEqAt] <- 1
    ttemp <- table(allAges, idtemp)
    temp <- c(ttemp[, 2])
    recTab[rownames(ttemp), at] <- temp
  }
  
  # Cumulative tables:
  cumTab <- recTab * 0
  for (at in ageTypes) {
    cumTab[, at] <- rev(cumsum(rev(recTab[, at])))
  }
  
  Nx <- cumTab[, "D"] + cumTab[, "C"] - cumTab[, "F"]
  Dx <- recTab[, "D"]
  
  idDead <- which(recTab[, "D"] > 0)
  nDead <- length(idDead)
  pleNB <- cumprod(1 - c(Dx / Nx)[idDead])
  
  # ==================== #
  # ==== BOOTSTRAP: ====
  # ==================== #
  bootTab <- matrix(0, nDead, nboot)
  for (iboot in 1:nboot) {
    idboot <- sample(1:n, size = n, replace = TRUE)
    ageLastb <- ageLast[idboot]
    ageFirstb <- ageFirst[idboot]
    departTypeb <- departType[idboot]
    
    # create identities and age list:
    if (is.null(ageFirst)) {
      allAgesb <- sort(ageLastb)
      allAgesIdb <- departTypeb
    } else {
      idAgeFirstb <- which(ageFirstb > min(ageLastb))
      allAgesb <- c(ageFirstb[idAgeFirstb], ageLastb)
      allAgesIdb <- c(rep("F", length(idAgeFirstb)), departTypeb)
    }
    
    unAllAgesb <- sort(unique(allAgesb))
    nagesb <- length(unAllAgesb)
    idmiss <- which(!unAllAges %in% unAllAgesb)
    allAgesb <- c(allAgesb, unAllAges[idmiss])
    allAgesIdb <- c(allAgesIdb, rep("N", length(idmiss)))
    
    ageTypesb <- unique(allAgesIdb)
    ntypesb <- length(ageTypesb)
    idsort <- sort.int(allAgesb, index.return = TRUE)$ix
    allAgesb <- allAgesb[idsort]
    allAgesIdb <- allAgesIdb[idsort]
    nAllAgesb <- length(allAgesb)
    
    # ageNames <- sprintf("a%s", gsub("[[:space:]]", "0", 
    #                                 format(1:nages, scientific = F)))
    
    # Count by type:
    recTab <- matrix(0, nages, NallTypes, 
                     dimnames = list(ageNames, allAgeTypes))
    for (at in ageTypes) {
      idtemp <- rep(0, nAllAges)
      idEqAt <- which(allAgesId == at)
      idtemp[idEqAt] <- 1
      ttemp <- table(allAges, idtemp)
      temp <- c(ttemp[, 2])
      recTab[rownames(ttemp), at] <- temp
    }
    
    # Cumulative tables:
    cumTab <- recTab * 0
    for (at in ageTypes) {
      cumTab[, at] <- rev(cumsum(rev(recTab[, at])))
    }
    
    Nx <- cumTab[, "D"] + cumTab[, "C"] - cumTab[, "F"]
    Dx <- recTab[, "D"]
    
    indDeadb <- rep(0, nages)
    indDeadb[which(recTab[, "D"] > 0)] <- 1
    pleb <- cumprod(((1 - c(Dx / Nx))^indDeadb)[idDead])
    bootTab[, iboot] <- pleb
  }
  
  # Calculate CIs:  # Calculate CIs:
  pleci <- t(apply(bootTab, 1, quantile, c(alpha / 2, 1 - alpha / 2),
                   na.rm = TRUE))
  colnames(pleci) <- c("Lower", "Upper")
  
  # Add age 0:
  Ages <- unAllAges[idDead]
  if (Ages[1] > 0) {
    Ages <- c(0, Ages)
    pleNB <- c(1, pleNB)
    pleci <- rbind(1, pleci)
  }
  
  pleboot <- data.frame(Ages = Ages, ple = pleNB, pleci)
  class(pleboot) <- "paramDemoPLECIs"
  return(pleboot)
}

# Plot Product-limit estimator:
plot.paramDemoPLE <- function(x, ...) {
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
  
  # Produce plots:
  op <- par(no.readonly = TRUE)
  
  par(mar = c(4, 4, 1, 1), mfrow = c(1, 1))
  plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
  
  lines(agev, x$ple, type = type, col = col, lwd = lwd)
  Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
  Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
  mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  mtext(xlab, side = 1, line = 2, cex = cex.lab)
  par(op)
}

# Plot Product-limit estimator CIs:
plot.paramDemoPLECIs <- function(x, ...) {
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
    if (length(args$col) < 2) {
      cols <- c(mean = "#800026", cis = "#FC4E2A")
    } else {
      cols <- args$col[1:2]
      names(cols) <- c("mean", "cis")
    }
  } else {
    cols <- c(mean = "#800026", cis = "#FC4E2A")
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
  
  # Produce plots:
  op <- par(no.readonly = TRUE)
  
  par(mar = c(4, 4, 1, 1), mfrow = c(1, 1))
  plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
  
  # polygon(c(agev, rev(agev)), c(x$Lower, rev(x$Upper)), col = cols["cis"], 
  #         border = NA)
  for (ici in c("Lower", "Upper")) {
    lines(agev, x[[ici]], type = type, col = cols["cis"], lwd = lwd, lty = 2)
  }
  lines(agev, x$ple, type = type, col = cols["mean"], lwd = lwd)
  Axis(xlim, side = 1, pos = ylim[1], cex.axis = cex.axis)
  Axis(ylim, side = 2, pos = xlim[1], las = 2, cex.axis = cex.axis)
  mtext(ylab, side = 2, line = mar[2] / 2, cex = cex.lab)
  mtext(xlab, side = 1, line = 2, cex = cex.lab)
  par(op)
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

# Simulate study:
SimulateStudy <- function(n, theta, model, shape, studyStart, studySpan,
                          birthStart, birthUncert = 0, dx = 1/365.25, 
                          provideSurv = TRUE, checkTheta = TRUE) {
  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)
  
  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }
  
  # Check studyStart date:
  if (is.character(studyStart)) {
    st <- try(as.Date(studyStart), silent = TRUE)
    stTest <- grepl("[[:alnum:]]{4}[[:punct:]]{1}[[:alnum:]]{2}[[:punct:]]{1}[[:alnum:]]{2}", studyStart)
    if (stTest) {
      stStop <- FALSE
      studyStart <- as.Date(studyStart)
    } else {
      stStop <- TRUE
    }
  } else if (class(studyStart) == "Date") {
    stStop <- FALSE
  } else {
    stStop <- TRUE
  }
  if (stStop) {
    stop("Argument 'studyStart' needs to be a date character string formatted as YYYY-MM-DD", call. = FALSE)
  }
  
  # Check birthStart date:
  if (is.character(birthStart)) {
    st <- try(as.Date(birthStart), silent = TRUE)
    stTest <- grepl("[[:alnum:]]{4}[[:punct:]]{1}[[:alnum:]]{2}[[:punct:]]{1}[[:alnum:]]{2}", birthStart)
    if (stTest) {
      stStop <- FALSE
      birthStart <- as.Date(birthStart)
    } else {
      stStop <- TRUE
    }
  } else if (class(birthStart) == "Date") {
    stStop <- FALSE
  } else {
    stStop <- TRUE
  }
  if (stStop) {
    stop("Argument 'birthStart' needs to be a date character string formatted as YYYY-MM-DD", call. = FALSE)
  }
  
  # Calculate parametric survival:
  if (provideSurv) {
    x <- seq(0, 1000, dx)
    Sx <- CalcSurv(theta, x, model = model, shape = shape, checkTheta = FALSE)
    idsx <- which(Sx >= 0.01)
    paramSurv <- list(calc = TRUE, 
                      surv = data.frame(x = x[idsx], surv = Sx[idsx]))
  } else {
    paramSurv <- list(calc = FALSE, surv = NA)
  }
  
  # Total number of sampled individuals:
  nAll <- n * 3
  
  # Set date of study end:
  studyEnd <- as.Date(julian(studyStart) + studySpan * 365.25, 
                      origin = "1970-01-01")
  
  # Simulate ages at death:
  agesDeath <- SampleRandAge(n = nAll, theta = theta, dx = dx, model = model, 
                             shape = shape)
  
  # Simulate times of birth:
  datesSeq <- seq(julian(birthStart), julian(studyEnd), dx)
  sampBirth <- sample(x = datesSeq, size = nAll, replace = TRUE)
  birthDates <- as.Date(sampBirth, origin = "1970-01-01")
  
  # Death dates:
  deathDates <- as.Date(sampBirth + agesDeath * 365.25, origin = "1970-01-01")
  
  # Find records within the study span:
  idincl <- which(deathDates >= studyStart)
  nincl <- length(idincl)
  idincl <- sort(idincl[sample(x = 1:nincl, size = n, replace = FALSE)])
  
  # subset data:
  agesDeath <- agesDeath[idincl]
  birthDates <- birthDates[idincl]
  deathDates <- deathDates[idincl]
  
  # Assign entry dates:
  entryDates <- birthDates
  entryDates[which(birthDates < studyStart)] <- studyStart
  
  # Find censored data:
  idCens <- which(deathDates > studyEnd)
  length(idCens)
  
  # Assign departure times:
  departDates <- deathDates
  departDates[idCens] <- studyEnd
  
  # Assign depart type:
  departType <- rep("D", n)
  departType[idCens] <- "C"
  
  # Find truncated records:
  idTrunc <- which(birthDates < studyStart)
  
  # Assign uncertain times of birth:
  minBirthDates <- maxBirthDates <- birthDates
  if (birthUncert > 0) {
    birthRange <- runif(n = length(idTrunc), min = 0, max = birthUncert * 365.25)
    minBirthDates[idTrunc] <- birthDates[idTrunc] - birthRange
    maxBirthDates[idTrunc] <- birthDates[idTrunc] + birthRange
  }
  maxBirthDates[idTrunc][which(maxBirthDates[idTrunc] > studyStart)] <- studyStart
  
  # Make data frame:
  output <- list(data = data.frame(ID = idincl, Birth.Date = birthDates, 
                                   Min.Birth.Date = minBirthDates,
                                   Max.Birth.Date = maxBirthDates, 
                                   Entry.Date = entryDates, 
                                   Depart.Date = departDates,
                                   Depart.Type = departType, 
                                   stringsAsFactors = FALSE),
                 surv = paramSurv)
  
  # Return output:
  return(output)
}

# ================================== END ===================================== #