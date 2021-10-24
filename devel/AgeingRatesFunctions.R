DefineAgeingRate <- function(theta, x, model = "GO", shape = "simple") {
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

# Test Weibull consistency:
CalcAR1 <- DefineAgeingRate(theta, x, model = "WE", shape = "simple") 

CalcAR2 <- DefineAgeingRate(theta, x, model = "WE", shape = "Makeham") 

CalcAR3 <- DefineAgeingRate(theta, x, model = "WE", shape = "bathtub")

ar1 <- CalcAR1(theta = c(b0 = 1.2, b1 = 0.5), x = 10)

ar2 <- CalcAR2(theta = c(c = 0, b0 = 1.2, b1 = 0.5), x = 10)

ar3 <- CalcAR3(theta = c(a0 = -Inf, a1 = 0, c = 0, b0 = 1.2, b1 = 0.5), x = 10)

cat(sprintf("simple:  %s\nMakeham: %s\nbathtub: %s\n", ar1, ar2, ar3))

# Test Gompertz consistency:
CalcAR1 <- DefineAgeingRate(theta, x, model = "GO", shape = "simple") 

CalcAR2 <- DefineAgeingRate(theta, x, model = "GO", shape = "Makeham") 

CalcAR3 <- DefineAgeingRate(theta, x, model = "GO", shape = "bathtub")

ar1 <- CalcAR1(theta = c(b0 = -1, b1 = 0.15), x = 10)

ar2 <- CalcAR2(theta = c(c = 0, b0 = -1, b1 = 0.15), x = 10)

ar3 <- CalcAR3(theta = c(a0 = -Inf, a1 = 0, c = 0, b0 = -1, b1 = 0.15), x = 10)

cat(sprintf("simple:  %s\nMakeham: %s\nbathtub: %s\n", ar1, ar2, ar3))

# Test Logistic consistency:
CalcAR1 <- DefineAgeingRate(theta, x, model = "LO", shape = "simple") 

CalcAR2 <- DefineAgeingRate(theta, x, model = "LO", shape = "Makeham") 

CalcAR3 <- DefineAgeingRate(theta, x, model = "LO", shape = "bathtub")

ar1 <- CalcAR1(theta = c(b0 = -1, b1 = 0.15, b2 = 0.5), x = 10)

ar2 <- CalcAR2(theta = c(c = 0, b0 = -1, b1 = 0.15, b2 = 0.5), x = 10)

ar3 <- CalcAR3(theta = c(a0 = -Inf, a1 = 0, c = 0, b0 = -1, b1 = 0.15, 
                         b2 = 0.5), x = 10)

cat(sprintf("simple:  %s\nMakeham: %s\nbathtub: %s\n", ar1, ar2, ar3))
