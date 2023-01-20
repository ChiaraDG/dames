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

### Complete Case Analyses

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

```
mod.wee    <- wee_mm(mean.formula = Y ~ X + time + Z + X:time,
                      t.formula = ~1,
                      id = id, weights = 1/samp.probs,
                      data = exampledat)
# print the results
summary(mod.wee)
#> Warning in summary.TwoPhaseMM(mod.wee): When performing a weighted likelihood analysis (by specifying the
#> samp.probi argument), robust standard errors are reported. Model based standard errors will not be correct and
#> should not be used.
#> 
#> Class:
#> TwoPhaseMM
#> 
#> Call:
#> wee_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1, 
#>     id = id, data = exampledat, weights = 1/samp.probs)
#> 
#> Information Criterion:
#>       AIC        BIC     logLik   Deviance  
#> 1845.1127  1864.6576  -916.5563  1833.1127  
#> 
#> Marginal Mean Model Parameters:
#>              Estimate Robust SE Chi Square  Pr(>Chi)
#> (Intercept) -2.157193  0.211541   103.9898 < 2.2e-16
#> X            0.904386  0.354652     6.5028  0.010770
#> time         0.189655  0.071099     7.1154  0.007642
#> Z            0.243668  0.249912     0.9507  0.329553
#> X:time       0.205454  0.143771     2.0421  0.152994
#> 
#> Dependence Model Parameters:
#>                   Estimate Robust SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  1.94548   0.21928     78.715 < 2.2e-16
#> 
#> Number of clusters:             192 
#> Maximum cluster size:           6 
#> Convergence status (nlm code):  1 
#> Number of iterations:           33
```

### Full-Cohort Analyses

- Multiple Imputation. Two MI methods are implemented: the indirect MI and the direct MI. Both methods support binary expensive covariates only. The indirect MI can be used only for NSA designs, whereas the direct MI can be used for both NSA and RDS designs.

```
mod.mi2 <- mi_mm(mean.formula = Y ~ X + time + Z + X:time,
                 t.formula = ~1,
                 lv.formula = NULL, id = "id", X = "X",
                 data = exampledat,
                 marg.exp.formula = X ~ Z, M = 5, method = "direct")
summary(mod.mi2)
#> 
#> Class:
#> TwoPhaseMM
#> 
#> Call:
#> mi_mm(mean.formula = Y ~ X + time + Z + X:time, lv.formula = NULL, 
#>     t.formula = ~1, id = "id", data = exampledat, sampled = "Sampled", 
#>     X = "X", M = 5, marg.exp.formula = X ~ Z, method = "direct")
#> 
#> Information Criterion:
#> [1]  NULL
#> 
#> Marginal Mean Model Parameters:
#>              Estimate  Model SE Chi Square  Pr(>Chi)
#> (Intercept) -2.217657  0.154648   205.6371 < 2.2e-16
#> X            1.124937  0.309988    13.1694 0.0002846
#> time         0.179099  0.052639    11.5765 0.0006679
#> Z            0.278590  0.177457     2.4646 0.1164373
#> X:time       0.169666  0.106359     2.5447 0.1106633
#> 
#> Dependence Model Parameters:
#>                   Estimate Model SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  2.10573  0.14999      197.1 < 2.2e-16
#> 
#> Number of clusters:             500 
#> Maximum cluster size:           6 
#> Convergence status (nlm code):  1 
#> Number of iterations:           NA
```

- Sieve Maximum Likelihood Estimator. The method can be used regardless of the type of expensive exposure (i.e., nominal or continuous) and/or study design.

```
# set up P(X|Z)
Bspline_Z <- cbind(as.numeric(exampledat$Z == 0), 
                   as.numeric(exampledat$Z == 1))
n_sieve             <- ncol(Bspline_Z)
colnames(Bspline_Z) <- paste("bs", 1:n_sieve, sep="")
exampledat          <- cbind(exampledat, Bspline_Z)
# estimate the parameters using SMLE
mod.smle <- smle_mm(mean.formula = Y ~ X + time + Z + X:time,
                    t.formula = ~1, lv.formula = NULL,
                    Y = "Y", X = "X", Z = "Z", Time = "time",
                    id = "id", Bspline = Bspline_Z, n_sieve = n_sieve, no_SE = FALSE,
                    data = exampledat)
summary(mod.smle)
#> 
#> Class:
#> TwoPhaseMM
#> 
#> Call:
#> smle_mm(Y = "Y", X = "X", Z = "Z", Time = "time", Bspline = Bspline_Z, 
#>     n_sieve = n_sieve, mean.formula = Y ~ X + time + Z + X:time, 
#>     t.formula = ~1, lv.formula = NULL, id = "id", data = exampledat, 
#>     no_SE = FALSE)
#> 
#> Information Criterion:
#>    logLik  
#> -932.8202  
#> 
#> Marginal Mean Model Parameters:
#>              Estimate  Model SE Chi Square  Pr(>Chi)
#> (Intercept) -2.192521  0.152959   205.4635 < 2.2e-16
#> X            0.978274  0.307756    10.1043  0.001479
#> time         0.160175  0.050374    10.1104  0.001474
#> Z            0.322854  0.179023     3.2523  0.071321
#> X:time       0.223221  0.116406     3.6772  0.055160
#> 
#> Dependence Model Parameters:
#>                   Estimate Model SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  2.10989  0.15045     196.67 < 2.2e-16
#> 
#> Number of clusters:             500 
#> Maximum cluster size:           5 
#> Convergence status (nlm code):  1 
#> Number of iterations:           18
```
