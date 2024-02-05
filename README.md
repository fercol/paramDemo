# paramDemo: Parametric and non-parametric demographic functions

`paramDemo` is an R package designed to explore age-specific survival, mortality, and fertility from parametric and non-parametric models.

## What's in paramDemo?

The package includes several functions to calculate age-specific survival, mortality, and fertility, as well as summary statistics from these demographic rates. It includes:  

- Functions to calculate age-specific survival and mortality from parametric mortality models following the conventions in package [BaSTA](https://github.com/fercol/BaSTA);
- Functions to calculate age-specific fertility from parametric models following the conventions in package [BaFTA](https://github.com/fercol/BaFTA);
- Functions to calculate non-parametric survival such as product limit estimators and life tables;
- Functions to calculate summary statistics such as remaining life expectancy, ageing rates, lifespan inequality and equality, age at maximum fecundity;
- Random simulation of longevities from age-specific mortality models;
- Functions to calculate and visualize stable population life history variables calculated from parametric models of age-specific survival and fertility.

## How to install paramDemo?
To install `paramDemo` from GitHub, type the following lines of code on the R console:

```R
# Install and load 'devtools':
install.packages("devtools")
library(devtools)

# Install paramDemo:
install_git("https://github.com/fercol/paramDemo", subdir = "pkg/")
```

