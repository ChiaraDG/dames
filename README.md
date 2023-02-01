# dames: Designs and Analysis Methods for Efficient Two-Phase Studies with a Binary Longitudinal Outcome

The `dames` package allows users to design and analyse cost and resource efficient epidemiological studies. The designs in the package can be used in scenarios where users have information on a correlated binary outcome and a set of confounders, and they are interested in a novel exposure that is expensive to collect. You can install the development version from GitHub with:

```
install.packages("devtools")
devtools::install_github("ChiaraDG/dames")
```

Once the package is installed, it can be loaded using:

```
library(dames)
```

## Designs

There are two classes of designs implemented in the `dames` package: the none-some-all design (NSA) and the residual dependent sampling (RDS) design. Below is an example of how the function `nsa_design()` can be used. Here we want to sample 85 people from the none stratum, 100 people from the some stratum and 15 people. Note, since within each of the three strata independent Bernoulli sampling is performed, it is not unusual to have a total number of sampled subjects being different than 200. Thus, in this README file we are going to set a seed for reproducibility.

```
set.seed(202311)
design <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
```

Let us assume that the expensive covariate was not observed in those patients who were not selected by `nsa_design()`. To do so, we set the values of X to missing if subject ID is not in the `sampled.id` section of `nsa_design()`:

```
# Assume X is not collected for subjects not sampled
exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)
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
     id = id, data = exampledat, samp.probs = design$sample.probs)
#>
#> Information Criterion:
#>      AIC        BIC     logLik   Deviance  
#> 285.6909   301.1392  -136.8454   273.6909  
#>
#> Marginal Mean Model Parameters:
#>              Estimate  Model SE Chi Square  Pr(>Chi)
#> (Intercept) -1.761520  0.278010    40.1470 2.356e-10
#> X            0.625321  0.489145     1.6343   0.20111
#> time         0.037619  0.102853     0.1338   0.71455
#> Z            0.684270  0.338074     4.0967   0.04297
#> X:time       0.270911  0.184323     2.1602   0.14163
#>
#> Dependence Model Parameters:
#>                   Estimate Model SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  1.87276  0.28595     42.891 5.787e-11
#>
#> Number of clusters:             97 
#> Maximum cluster size:           6 
#> Convergence status (nlm code):  1 
#> Number of iterations:           29
```
  
- Weighted Maximum Likelihood Estimator. The method can be used regardless of the type of expensive exposure (i.e., nominal or continuous). The method can be implemented as long as every subject has a non-zero probability of being sampled; thus, it cannot be implemented for the RDS designs current implemented in `dames`. The argument `weights` in the `wee_mm()` call needs to be a vector of length equal to the number of rows in the data.frame `data`

```
# create the vector of sampling probabilities for each subject
samp.probs <- ifelse(design$membership == "None", design$sample.probs[1], 
                     ifelse(design$membership == "Some", design$sample.probs[2], 
                            design$sample.probs[3]))
# number of times each subject is observed
time_i     <- exampledat %>% group_by(id) %>% summarise(mi = n()) %>% pull()
# repeat the sampling probability for each subject based on the number
# of times a subject is observed
exampledat$samp.probs <- rep(samp.probs, time_i)
samp.probs            <- exampledat[!is.na(exampledat$X),"samp.probs"]

mod.wee    <- wee_mm(mean.formula = Y ~ X + time + Z + X:time,
                      t.formula = ~1,
                      id = id, weights = 1/samp.probs,
                      data = exampledat)
# print the results
summary(mod.wee)
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
#>  945.5626   961.0109  -466.7813   933.5626  
#>
#> Marginal Mean Model Parameters:
#>              Estimate Robust SE Chi Square  Pr(>Chi)
#> (Intercept) -1.892993  0.300329    39.7286 2.918e-10
#> X            0.636161  0.514212     1.5306    0.2160
#> time         0.050214  0.102960     0.2379    0.6258
#> Z            0.812207  0.387953     4.3830    0.0363
#> X:time       0.262024  0.181104     2.0933    0.1479
#>
#> Dependence Model Parameters:
#>                   Estimate Robust SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  1.75827   0.30514     33.203 8.301e-09
#>
#> Number of clusters:             97 
#> Maximum cluster size:           6 
#> Convergence status (nlm code):  1 
#> Number of iterations:           31
#> Warning message:
#> In summary.TwoPhaseMM(mod.wee) :
#>   When performing a weighted likelihood analysis (by specifying the weights argument), robust standard #> errors are reported. Model based standard errors will not be correct and should not be used.
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
#>     t.formula = ~1, id = "id", data = exampledat, X = "X", M = 5, 
#>     marg.exp.formula = X ~ Z, method = "direct")
#>
#> Information Criterion:
#> [1]  NULL
#>
#> Marginal Mean Model Parameters:
#>              Estimate  Model SE Chi Square Pr(>Chi)
#> (Intercept) -2.112184  0.150236   197.6598  < 2e-16
#> X            0.967794  0.513037     3.5585  0.05924
#> time         0.165287  0.063028     6.8772  0.00873
#> Z            0.278927  0.225904     1.5245  0.21694
#> X:time       0.228862  0.196857     1.3516  0.24500
#>
#> Dependence Model Parameters:
#>                   Estimate Model SE Chi Square  Pr(>Chi)
#> gamma:(Intercept)  2.12042  0.14773     206.02 < 2.2e-16
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
