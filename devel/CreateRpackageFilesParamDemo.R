# Source files:
# source("~/FERNANDO/PROJECTS/4.PACKAGES/CreateRproj/pkg/R/CreateRproj.R")
library(CreateRproj)
# Main directory:
mainDir <- "~/FERNANDO/PROJECTS/4.PACKAGES/"

# Package name:
pkgName <- "paramDemo" 

# Code name:
scriptFile <- "paramDemo"

# Create .Rd (help) files (requires a code file):
CreateRdFiles(pkgName = pkgName, mainDir = mainDir, scriptFile = scriptFile,
              authorNames = "Fernando Colchero", 
              authorEmails = "colchero@imada.sdu.dk")

# Create namespace:
CreateNamespace(pkgName = pkgName, mainDir = mainDir, scriptFile = scriptFile)

# Create description file:
CreatePkgDescrip(pkgName = pkgName, mainDir = mainDir, title = "Parametric and non-parametric demographic functions and applications.", version = "1.0.0", authors = "Fernando Colchero", maintainer = "Fernando Colchero <colchero@imada.sdu.dk>")