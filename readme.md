# hyb_niches

## goal

## dependencies 

### VS Code
To run in VS Code, make sure you have the extensions "Julia", "R", and "Quarto". 

### julia

to instantiate / load packages from `\hyb_niches` directory:    

`julia --project`

```julia
using Pkg
# install the packages listed in the environment
Pkg.instantiate()
```

### R

packages:
```R
library(devtools)
list_of_packages <- c(
  # Spatial / Mapping
  "sf", "ggspatial",
  # Visualization
  "ggplot2", "viridis", "gridExtra",
  # Data manipulation
  "dplyr", "tidyr", "stringr",
  # Biodiversity data
  "ridigbio", "gatoRs"
)

new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if (length(new.packages)) {
  cat("\nInstalling missing CRAN packages:\n")
  print(new.packages)
  install.packages(new.packages)
} else {
  cat("\nAll required CRAN packages already installed.\n")
}

## ---- Load All CRAN Packages ----

# Try loading all required CRAN packages.
loaded <- sapply(list_of_packages, require, character.only = TRUE)

if (any(!loaded)) {
  cat("\nWARNING: These CRAN packages failed to load:\n")
  print(list_of_packages[!loaded])
} else {
  cat("\nAll CRAN packages loaded successfully.\n")
}
```
