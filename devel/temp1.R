# Main function:
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
  class(outList) <- c("paramDemo", ifelse(SURV & FERT, "demoBoth", 
                           ifelse(SURV, "demoSurv", "demoFert")))
  return(outList)
}


dem1 <- CalcDemo(theta = c(-3, 0.1), beta = c(0.5, 0.02, 3), ageMatur = 2)
dem2 <- CalcDemo(beta = c(0.5, 0.02, 3), ageMatur = 2, type = 'fertility', maxAge = 100)
