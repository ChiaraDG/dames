# dames: Designs and Analysis Methods for Efficient Two-Phase Studies with a Binary Longitudinal Outcome

The `dames` package allows users to design and analyse cost and resource efficient epidemiological studies. The designs in the package can be used in scenarios where users have information on a correlated binary outcome and a set of confounders, and they are interested in a novel exposure that is expensive to collect. You can install the development version from GitHub with:

```
install.packages("devtools")
devtools::install_github("ChiaraDG/dames")
```

## Designs

There are two classes of designs implemented in the `dames` package: the none-some-all design (NSA) and the residual dependent sampling (RDS) design. Below is an example of how the function `nsa_design()` can be used. Here we want to sample 85 people from the none stratum, 100 people from the some stratum and 15 people. Note, since within each of the three strata independent Bernoulli sampling is performed, it is not unusual to have a total number of sampled subjects being different than 200.

```
design <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
```

## Analysis Methods

Inference methods can be divided into two groups: methods that only include subjects with complete data on outcome, expensive exposure and inexpensive covariates (complete-case analyses), and methods that include all subjects regardless of whether they have information on the expensive exposure (full-cohort analyses).

The `dames` package allows users to perform both complete-case and full-cohort analyses:
* Complete Case Analyses

  - Ascertainment Corrected Maximum Likelihood Estimator: the method can be used regardless of the type of expensive exposure (i.e., nominal or continuous). The method can be used for NSA designs only.
  
```
mod.acml    <- acml_mm(mean.formula = Y ~ X + time + Z + X:time,
                       t.formula = ~1,
                       id = id, samp.probs = design$sample.probs,
                       data = exampledat)

# print the results
summary(mod.acml)
#> 
#> Class:
#> TwoPhaseMM
#> 
#> Call:
#> acml_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1, 
#>     id = id, data = exampledat, samp.probs = design$sample.probs)
#> 
#> Information Criterion:
#>       AIC        BIC     logLik   Deviance  
#>  515.4883   535.0333  -251.7442   503.4883  
#> 
#> Marginal Mean Model Parameters:
#>              Estimate  Model SE Chi Square Pr(>Chi)
#> (Intercept) -2.100535  0.210039   100.0143  < 2e-16
#> X            0.804936  0.350429     5.2762  0.02162
#> time         0.179964  0.072248     6.2047  0.01274
#> Z            0.267792  0.234088     1.3087  0.25263
#> X:time       0.233289  0.135775     2.9522  0.08576
#> 
#> Dependence Model Parameters:
#>                   Estimate Model SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  1.98256  0.20576      92.84 < 2.2e-16
#> 
#> Number of clusters:             192 
#> Maximum cluster size:           6 
#> Convergence status (nlm code):  1 
#> Number of iterations:           31
```
  
  - Weighted Maximum Likelihood Estimator. The method can be used regardless of the type of expensive exposure (i.e., nominal or continuous). The method can be implemented as long as every subject has a non-zero probability of being sampled; thus, it cannot be implemented for the RDS designs described in this R package.

* Full-Cohort Analyses

  - Multiple Imputation. Two MI methods are implemented: the indirect MI and the direct MI. Both methods support binary expensive covariates only. The indirect MI can be used only for NSA designs, whereas the direct MI can be used for both NSA and RDS designs.

  - Sieve Maximum Likelihood Estimator. The method can be used regardless of the type of expensive exposure (i.e., nominal or continuous) and/or study design.


