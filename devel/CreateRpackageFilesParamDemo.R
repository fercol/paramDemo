# Source files:
# source("~/FERNANDO/PROJECTS/4.PACKAGES/CreateRproj/pkg/R/CreateRproj.R")
library(CreateRproj)
# Main directory:
mainDir <- "~/FERNANDO/PROJECTS/4.PACKAGES/"

# Package name:
pkgName <- "paramDemo" 

# Code name:
codeFile <- "paramDemo"

# Create .Rd (help) files (requires a code file):
CreateRdFiles(pkgName = pkgName, mainDir = mainDir, codeFile = codeFile,
              authorNames = "Fernando Colchero", 
              authorEmails = "colchero@imada.sdu.dk")

# Create namespace:
CreateNamespace(pkgName = pkgName, mainDir = mainDir, codeFile = codeFile)

