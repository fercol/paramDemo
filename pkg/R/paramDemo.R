# ============================ CODE METADATA ================================= #
# PACKAGE: paramDemo
# AUTHOR: Fernando Colchero
# DATE CREATED: 2020-06-01
# DATE MODIFIED: 2021-10-14
# DESCRIPTION: Functions to extract demographic information from parametric
#              mortality and Fertility models, summary statistics (e.g.
#              ageing rates, life expectancy, lifespan equality, etc.)
#              and life table calculation from census data.
# COMMENTS: Created a function to calculate CIs of life table elements and
#           plotting functions for life tables and life table CIs.
# ============================ START CODE ==================================== #
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
  if (!modelFert %in% sprintf("M%s", 1:3)) {
    stop("Wrong 'modelFert' specification.\n",
         "       Alternatives are 'M1', 'M2', or 'M3'.", call. = FALSE)
    
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
      lowTh[3] == -Inf
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
    stop(sprintf("Some theta parameters are below their lower bound.\n %s.\n",paste(sprintf("min(%s) = %s", nameTh, lowTh), collapse = ", ")),
         call. = FALSE)
    
  }
  
  defaultTheta  <- list(theta = theta, p = nTh, name = nameTh,
                        low = lowTh)
  attr(defaultTheta, "model") = model
  attr(defaultTheta, "shape") = shape
  return(defaultTheta)
}

# Fertility:
.SetBeta <- function(beta, modelFert = "M1") {
  if (is.null(beta)) {
    stop("Missing 'beta' parameter vector or matrix.\n", call. = FALSE)
  } 
  if (modelFert == "M1") {
    nBe <- 3
    lowBe <- c(0, 0, 0)
  } else if (modelFert == "M2") {
    nBe <- 4
    lowBe <- c(0, 0, -Inf, 0)
  } else if (modelFert == "M3") {
    nBe <- 4
    lowBe <- c(0, 0, 0, -Inf)
  }
  nameBe <- sprintf("b%s", 1:nBe - 1)
  if (is.matrix(beta)) {
    nbeUser <- ncol(beta)
    clBe <- "matrix"
    stBe <- "columns"
  } else if (is.numeric(beta)) {
    nbeUser <- length(beta)
    clBe <- "vector"
    stBe <- "elements"
  } else {
    stop("Parameters beta should either be of class matrix\nor a numeric vector.\n", call. = FALSE)
  }
  if (nbeUser != nBe) {
    stop(sprintf("The beta %s should have %s %s.\n", clBe, nBe, stBe),
         call. = FALSE)
  } else {
    if (is.matrix(beta)) {
      colnames(beta) <- nameBe
    } else {
      names(beta) <- nameBe
    }
  }
  # check if beta parameters conform to their support:
  if (is.matrix(beta)) {
    BETLOW <- all(sapply(1:nBe, function(bi) {
      bl <- all(beta[, bi] >= lowBe[bi])
    }))
  } else {
    BETLOW <- all(beta >= lowBe)
  }
  if(!BETLOW) {
    stop(sprintf("Some beta parameters are below their lower bound.\n %s.\n",paste(sprintf("min(%s) = %s", nameBe, lowBe), collapse = ", ")),
         call. = FALSE)
    
  }
  defaultBeta  <- list(beta = beta, p = nBe, name = nameBe,
                        low = lowBe)
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
          x^(theta[, "b0"] - 1)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + theta[, "b0"] * theta[, "b1"]^theta[, "b0"] *
          x^(theta[, "b0"] - 1)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] +
          theta[, "b0"] * theta[, "b1"]^theta[, "b0"] *
          x^(theta[, "b0"] - 1)
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
          x^(theta["b0"] - 1)
      }
    } else if (shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + theta["b0"] * theta["b1"]^theta["b0"] *
          x^(theta["b0"] - 1)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] +
          theta["b0"] * theta["b1"]^theta["b0"] *
          x^(theta["b0"] - 1)
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
.DefineFertilityNumeric <- function(modelFert = "M1") {
  if (modelFert == "M1") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2)
      return(fert)
    }
  } else if (modelFert == "M2") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2) /
        (1 + exp(-(beta["b3"] * (x - beta["b2"]))))
      return(fert)
    }
  } else if (modelFert == "M3") {
    fertfun <- function(beta, x) {
      fert <- beta["b0"] * exp(-beta["b1"] * (x - beta["b2"])^2 +
                   beta["b3"] * 1/(x + 1))
      return(fert)
    }
  }
  return(fertfun)
}

.DefineFertilityMatrix <- function(modelFert = "M1") {
  if (modelFert == "M1") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2)
      return(fert)
    }
  } else if (modelFert == "M2") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2) /
        (1 + exp(-(beta[, "b3"] * (x - beta[, "b2"]))))
      return(fert)
    }
  } else if (modelFert == "M3") {
    fertfun <- function(beta, x) {
      fert <- beta[, "b0"] * exp(-beta[, "b1"] * (x - beta[, "b2"])^2 +
                   beta[, "b3"] * 1/(x + 1))
      return(fert)
    }
  }
  return(fertfun)
}

# ------------- #
# AGEING RATES:
# ------------- #
.DefineAgeingRate <- function(model = "GO", shape = "simple") {
  # -------------------------------------
  # Exponential (i.e. constat mortality):
  # -------------------------------------
  if (model == "EX") {
    # ------------
    # Exponential:
    # ------------
    ageingRate <- function(theta, x) rep(0, length(x))
    # ---------
    # Gompertz:
    # ---------
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
    # --------
    # Weibull:
    # --------
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
    # ---------
    # Logistic:
    # ---------
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
          (theta["b1"] * (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] *
                            (exp(theta["b1"] * x) - 1)) -
             theta["b2"] * exp(theta["b0"] + theta["b1"] * x)) /
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


# ------------------------- #
# PACE AND SHAPE FUNCTIONS:
# ------------------------- #
# Life expectancy, lifespan inequality, lifespan equality, Gini, coeff. of
# variation:
.CalcPaceShape <- function(theta, x, dx, survFunc) {
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

# --------------- #
# USER FUNCTIONS:
# --------------- #
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
Calcfert <- function(beta, x, modelFert = "M1", checkBeta = TRUE) {

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

# Main demographic function:
CalcDemo <- function(theta = NULL, beta = NULL, x = NULL, dx = NULL, 
                     model = "GO", shape = "simple", modelFert = "M1", 
                     type = "both", minSx = 0.01, summarStats = TRUE, 
                     ageMatur = 0, maxAge = NULL, agesAR = NULL, 
                     SxValsAR = NULL) {
  
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
      ageingRates <- cbind(CalcAgeingRate(theta = theta, x = agesAR, 
                                          model = model, shape = shape), 
                           Surv = round(SxValsAR, 4))
      
      # Pace-shape:
      paceShape <- .CalcPaceShape(theta, xInt, dxInt, survFunc = .CalcSurv)
      
      # Logical for calculation of summary statistics:
      calc <- TRUE
    } else {
      ageingRates <- NA
      paceShape <- NA
      calc <- FALSE
    }
    
    # List of outputs:
    survList <- list(functs = data.frame(age = x, mort = mort, surv = surv, 
                                         pdf = dens, cumhaz = cumhaz), 
                     summStats = list(calculated = calc, 
                                      ageingRates = ageingRates,
                                      paceShape = paceShape), 
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
        stop("Missing argument 'x' for vector of ages, or argument for maximum age 'maxAge'. Provide either.", call. = FALSE)
      } else {
        if (is.null(dx)) dx <- 0.01
        x <- seq(ageMatur, maxAge, dx)
      }
    }
    
    # Calculate age-specific fertility:
    fert <- .CalcFert(beta, x[which(x >= ageMatur)] - ageMatur)
    
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
                     sumStats = list(calculated = calc, ageMaxFert),
                     settings = list(beta = beta, modelFert = modelFert), 
                     analyzed = TRUE)
  } else {
    fertList <- list(analyzed = FALSE)
  }
  
  # ========================== #
  # III) CREATE OUTPUT OBJECT: 
  # ========================== #
  outList <- list(surv = survList, fert = fertList)
  class(outList) <- ifelse(SURV & FERT, "demoBoth", 
                           ifelse(SURV, "demoSurv", "demoFert"))
  return(outList)
}

# Function to calculate remaining life expectancy:
CalcRemainLifeExp <- function(theta, x = NULL, dx = NULL, xmax = NULL,
                              atAllAges = FALSE, model = "GO", shape = "simple",
                              checkTheta = TRUE) {
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

# Ageing rates:
CalcAgeingRate <- function(theta, x, model = "GO", shape = "simple",
                           checkTheta = TRUE) {
  # Verify model and shape:
  .VerifySurvMod(model = model, shape = shape)
  
  # Extract theta attributes:
  if (checkTheta) {
    thetaAttr <- .SetTheta(theta, model = model, shape = shape)
    theta <- thetaAttr$theta
  }
  
  # Ageing rate function:
  arFun <- .DefineAgeingRate(model = model, shape = shape)
  
  # Calculate ageing rate:
  ar <- cbind(Age = x, AR = arFun(theta = theta, x = x))
  rownames(ar) <- NULL
  
  return(ar)
}

# Age-specific survival probability:
CalcSurvProbs <- function(demo, dx = 1) {
  if (class(demo) != "demoSurv") {
    stop("Object 'demo' should be of class 'demoSurv'.\nCreate object demo with function CalcDemo().", call. = FALSE)
  }
  if (length(demo$age) == 1) {
    stop("Age-specific probabilities cannot be calculated from a single age.\nIncrease the age vector on CalcDemo().", call. = FALSE)

  }
  x <- demo$age
  mux <- demo$mort
  cumh <- demo$cumhaz
  Sx <- demo$surv
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

# Calculate age at maximum Fertility:
CalcAgeMaxFert <- function(beta, modelFert = "M1", ageMatur = 0, 
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
  if (modelFert == "M1") {
    xm <- beta["b2"]
    dd <- 0
    ii <- 0
  } else if (modelFert %in% c("M2", "M3")) {
    if (modelFert == "M2") {
      dfdx <- function(x, beta) {
        w <- x - beta["b2"]
        u <- exp(-beta["b1"] * w^2)
        z <- exp(-beta["b3"] * w)
        v <- 1 + z
        du <- -2 * beta["b1"] * w * u
        dv <- -beta["b3"] * z
        df <- (du * v - u * dv) / v^2
        return(df)
      }
      dfdx2 <- function(x, beta) {
        w <- x - beta["b2"]
        u <- exp(-beta["b1"] * w^2)
        z <- exp(-beta["b3"] * w)
        v <- 1 + z
        du <- -2 * beta["b1"] * w * u
        dv <- -beta["b3"] * z
        ddu <- (-2 * beta["b1"] + (-2 * beta["b1"] * w)^2) * u
        ddv <- beta["b3"]^2 * z
        df <- (du * v - u * dv) / v^2
        df2 <- (ddu * v - u * ddv) / v^2 - 2 * dv / v * df
        return(df2)
      }
    } else {
      dfdx <- function(x0, beta) {
        x0^3 + x0^2 * (2 - beta["b2"]) + x0 * (1 - 2 * beta["b2"]) +
          beta["b3"] / (2 * beta["b1"]) - beta["b2"]
      }
      dfdx2 <- function(x0, beta) {
        3 * x0^2 + 2 * x0 * (2 - beta["b2"]) + 1 - 2 * beta["b2"]
      }
    }
    id0 <- which(sign(dfdx(xv[-length(xv)], beta)) != sign(dfdx(xv[-1], beta)))
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
  }
  maxFert <- .CalcFert(beta, xm)
  maxFertv <- c(xm + ageMatur, maxFert, dd, ii, ageMatur)
  names(maxFertv) <- c("Age", "maxFert", "error", "iterations", "ageMatur")
  return(maxFertv)
}

# Sampling random ages at death:
SampleRandAge <- function(n, theta, dx = 0.001, model = "GO", shape = "simple",
                          minSx = 0.0001) {

  # Extract demographic functions:
  demo <- CalcDemo(theta, dx = dx, model = model, shape = shape, minSx = minSx,
                   summarStats = FALSE)

  # Extract CDF of ages at death:
  Fx <- 1 - demo$surv

  # Draw random uniform values:
  u <- runif(n)

  # Extract ages for Fx = u:
  uages <- demo$age[findInterval(u, Fx, rightmost.closed = TRUE)]

  # return random ages:
  return(uages)
}

# survival and average adult survival:
FindSilerPars <- function(theta, palpha) {
  theta[c(2, 3, 5)] <- abs(theta[c(2, 3, 5)])

  # Extract mortality, survival, etc:
  demotest <- CalcDemo(theta = theta, shape = "bathtub", minSx = 0.00001)

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

# Recursive optimization function:
myoptim <- function(par, fn, ...) {
  opt1 <- optim(par = par, fn = fn, ...)
  ntry <- 0
  while(opt1$convergence != 0 & ntry < 50) {
    ntry <- ntry + 1
    opt1 <- optim(par = opt1$par, fn = fn, ...)
  }
  return(opt1)
}

# Extract average juvenile and adult vital rates:
CalcAveDemo <- function(demo) {

  # Find life table probabilities
  probs <- CalcProbs(demo)

  # Find age at minimum mortality:
  alpha <- floor(demo$age[which(demo$mort == min(demo$mort))])

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

# Calculate life table:
CalcLifeTable <- function(ageLast, ageFirst = NULL, departType) {
  # Number of records:
  n <- length(ageLast)

  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, n)
  }

  # Unit age vector for that sex:
  agev <- 0:ceiling(max(ageLast))
  nage <- length(agev)

  # Outputs:
  Nx <- Dx <- ax <- rep(0, nage)
  for (xx in 1:nage) {
    # A) EXPOSURES:
    # Find how many entered the interval (including truncated):
    idNx <- which(ageFirst < agev[xx] + 1 & ageLast >= agev[xx])

    # Extract ages and departType:
    xf <- ageFirst[idNx]
    xl <- ageLast[idNx]
    dt <- departType[idNx]

    # proportion of truncation in interval:
    trp <- xf - agev[xx]
    trp[trp < 0] <- 0

    # proportion of censoring:
    cep <- agev[xx] + 1 - xl
    cep[cep < 0] <- 0
    cep[dt == "D"] <- 0

    # Calculate exposures:
    nexp <- 1 - trp - cep
    Nx[xx] <- sum(nexp)

    # B) DEATHS:
    # Calculate total deaths in the interval:
    idDx <- which(dt == "D" & xl < agev[xx] + 1)
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
  Tx <- rev(cumsum(rev(Lx)))

  # Remaining life expectancy after age x:
  ex <- Tx / lx
  ex[which(is.na(ex))] <- 0
  # (Note: follows on the correction for Lx)
  # ex <- Tx / Nx

  # Life-table:
  lt <- cbind(Ages = agev, Nx = Nx, Dx = Dx, lx = lx, px = px,
              qx = qx, Lx = Lx, Tx = Tx, ex = ex)
  class(lt) <- c("lifetable")
  return(lt)
}

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

# Plot lifetable:
plot.lifetable <- function(x, demorate = "lx", ...) {
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
plot.lifeTableCIs <- function(x, demorate = "lx", ...) {
  
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
