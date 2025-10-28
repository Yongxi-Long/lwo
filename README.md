
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lwo

<!-- badges: start -->
<!-- badges: end -->

The goal of lwo is to estimate Longitudinal Win Odds (lwo) from
longitudinal extension of the Probabilistic Index Model (PIM)

## Installation

You can install the development version of lwo from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Yongxi-Long/lwo")
```

## Example

The key function is *lwo()*, which fits a longitudinal probabilistic
index model and estimate longitudinal win odds for ordinal repeated
measurements.

We use example data set *dat_SID* from this package for illustration. It
is a simulated data set derived from the SID-GBS trial (Second
intravenous immunoglobulin dose in patients with Guillain-Barré syndrome
with poor prognosis (SID-GBS): a double-blind, randomised,
placebo-controlled trial Walgaard, ChristaWalgaard, Christa et al. The
Lancet Neurology). *dat_SID* contains both the wide format data.frame
and the long format data.frame, of the following variables:

- patient_ID: unique patient identifier
- age: patient age at entry
- pre_diarrhea: whether had preceding diarrhea before entry
- treat_ITT: treatment assignment under intention-to-treat
- GBS_DS: Guillain-Barré syndrome disability scale at week 0, 1, 2, 4,
  8, 12, 26

### The *lwo()* function

``` r
library(lwo)
# load the data set
data("dat_SID")
# extract the long format
dat_SID_long <- dat_SID$long_format
# fit a longitudinal PIM model
# we estimate the treatment effect quantified as win odds,
# adjusted for age, preceding diarrhea and baseline GBS-DS score
# we also model treatment-by-time interaction via a spline term
mod.lwo.spline <- lwo(
  GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline+
    treat_ITT:splines::ns(week,knots = 2),
  data = dat_SID_long,
  id.var = "patient_ID",
  visit.var = "week",
  time.vars = "week",
  corstr = "ar1",
  larger = FALSE
)
summary(mod.lwo.spline)
#> 
#> Call:
#> lwo(formula = GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline + 
#>     treat_ITT:splines::ns(week, knots = 2), data = dat_SID_long, 
#>     id.var = "patient_ID", visit.var = "week", time.vars = "week", 
#>     larger = FALSE, corstr = "ar1")
#> 
#> Coefficients:
#>                                          Estimate   Std.err   Wald Pr(>|W|)    
#> treat_ITT                               -0.206172  0.192056  1.152 0.283048    
#> age                                      0.008956  0.008613  1.081 0.298407    
#> pre_diarrhea                             0.194412  0.137678  1.994 0.157926    
#> GBS_DS_baseline                         -0.219022  0.052684 17.283 3.22e-05 ***
#> treat_ITT:splines::ns(week, knots = 2)1  2.016802  0.562635 12.849 0.000338 ***
#> treat_ITT:splines::ns(week, knots = 2)2  1.102107  0.250728 19.322 1.10e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Temporal Correlation Structure: ar1

newdata <- dat_SID_long |>
  dplyr::filter(patient_ID %in% c(1,4))
predict(object = mod.lwo.spline,
        newdata = newdata,
        id.var = "patient_ID",
        visit.var = "week",
        time.vars = "week",
        type = "link",
        conf.int = TRUE) |>
  exp() |>
  data.frame()
#>    lower.CI  estimate  upper.CI
#> 1 1.2862759 2.0277075 3.1965130
#> 2 1.1774330 1.7710568 2.6639668
#> 3 0.9537879 1.3855316 2.0127093
#> 4 0.6117448 0.9397939 1.4437600
#> 5 0.4413670 0.7156131 1.1602636
#> 6 0.2751556 0.4729790 0.8130276
```

### The *gen_data()* function

The package also provides a function, the *gen_data()*, to simulate
ordinal longitudinal data based on user inputs. The function internally
relies on the *genOrdCat()* function from the *simstudy* package.

``` r
# define the covariate distribution
# if covariate_def is left NULL, the function will use build-in distributions for (continuous) age variable and (binary) preceding diarrhea variable.
def <- simstudy::defData(varname = "sex", dist = "binary", formula = 0.5)
def <- simstudy::defData(def, varname = "age50", dist = "normal", formula = 0, variance = 100)
# define baseline category probabilities when all covariates are zero/reference values
baseprobs <- c(0.3,0.4,0.3)
# correlation matrix, must match the dimension of the visits
corMatrix <- gen_corMatrix(n_visits = 4,rho=0.6,corstr = "ar1")
dat <- gen_data(
   N = 200,
   baseprobs = baseprobs,
   covs_effects = c(sex = 0.5, age50 = -0.02),
   time_effects = c(0.1, 0.3, 0.5),
   trt_ratio = 1,
   time_trt_effects = c(0.05, 0.1, 0.15),
   visits = paste0("v", 1:4),
   corMatrix = corMatrix,
   covariate_def = def
 )
```

### The *calculate_win_odds()* function

The longitudinal ordinal data is generated from a proportional odds
model, in which the effects are specified as log odds ratio. As a
results, it is not straightforward to see what the true win odds is for
each visit after specifying the data generation paramters. This package
provides a function, the *calculate_win_odds()*, that use simulation to
compute the true (population-level) win odds at each visit, marginalized
with respect to the covariate distribution.

The input of this function is similar to that of the *gen_data()*
function, except that we have *N_approx* instead of *N* to specify the
number of iterations of the simulation to approximate the true win odds.

``` r
corMatrix <- gen_corMatrix(n_visits = 3,rho=0.6,corstr = "ar1")
estimands <- calculate_win_odds(N_approx = 1e4,
                               baseprobs = rev(c(0.06,0.11,0.12,0.50,0.21)),
                               covs_effects = c("age"=-0.005,
                                                "pre_diarrhea"= 0.23),
                               time_effects = c(0.6,1.2),
                               trt_ratio = 1,
                               time_trt_effects = c(0.4,0.8),
                               visits = visits <- c("week0","week4","week8"),
                               corMatrix = corMatrix)
#> Covariate distribution not supplied!
#> Will use build-in distributions for (continuous) age variable and (binary) preceding diarrhea variable.
estimands
#>    week0    week4    week8 
#> 1.000000 1.268990 1.649854
```
