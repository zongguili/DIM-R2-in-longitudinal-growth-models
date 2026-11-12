# DIM R Functions: Illustration Examples

This repository provides illustrative examples for applying the DIM R functions developed by Li and Ryu (2026) to compute R² measures for selected fixed and random effects in longitudinal growth and multilevel models. The accompanying paper uses empirical data from existing longitudinal study databases, whereas this repository uses a simulated longitudinal dataset for the illustrations.

## Repository files

- `r2lmer.R`: computes R² measures for models estimated with `lme4::lmer()`.
- `r2lme.R`: computes R² measures for models estimated with `nlme::lme()`.
- `r2lmeRmat.R`: computes R² measures for `lme()` models with more complex level-1 residual covariance structures.
- `bp_example.csv`: simulated longitudinal blood-pressure dataset used in the illustrations.
- `Illustration_examples.Rmd`: R Markdown file containing the model specifications, function calls, output, and interpretation examples.
- `Illustration_examples.pdf`: PDF version of the illustration for users who want to view the complete examples without running the R Markdown file.

## Requirements

The examples require the following R packages:

```r
install.packages(c("lme4", "nlme"))
```

Load the packages and DIM functions before running the examples:

```r
library(lme4)
library(nlme)

source("r2lmer.R")
source("r2lme.R")
source("r2lmeRmat.R")
```

## Reproducing the illustrations

1. Download this repository.
2. Keep `bp_example.csv`, the three DIM function files, and `Illustration_examples.Rmd` in the same working directory.
3. Open `Illustration_examples.Rmd` in RStudio.
4. Run the code chunks sequentially or knit the file to HTML.

If you only want to view the illustration and output, open `Illustration_examples.pdf`.

## Functions

The three functions use a consistent input structure:

- `model`: a fitted `lmer()` or `lme()` model object.
- `effectf`: the target fixed-effect term(s). Use `NULL` if no fixed effects are selected.
- `effectr`: the target random-effect term(s). Use `NULL` if no random effects are selected; use `"(Intercept)"` to select the random intercept.

See `Illustration_examples.Rmd` or `Illustration_examples.pdf` for complete examples.

## Reference

Li, Z., & Ryu, E. (2026). *R functions to compute effect size measures in longitudinal growth models*.[Manuscript under review].

