plot.paramDemo <- function(x, demofun = "all", ...) {
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
}

plot(dem1, lwd = 4, demofun = 'all')
plot(dem2, lwd = 4, xlim = c(0, 20))
