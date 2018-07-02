Introduction
------------

This document is a supplement to "Evaluating generalised additive mixed modelling strategies for dynamic speech analysis," relating specifically to the contents of Table 3 in Section 3.2. It presents code that illustrates (i) how the simulated data were generated and (ii) the models whose performance is summarised Table 3.

Preliminaries
-------------

The code below loads the relevant libraries.

``` r
library(ggplot2)
library(mgcv)
library(itsadug)
```

Data generation
---------------

The code in this section can be used to create data for either type I or type II simulations. Set the value of *type* to 1 for type I simulations and to 2 for type II simulations.

``` r
type = 2
```

The data for this set of simulations consist of simulated f2 trajectories loosely modelled after the diphthong /eI/. In this simulation, there are 50 different word types, each of which is represented by 10 trajectories. For type I simulations, each word is randomly assigned to one of two groups (A and B), but there is no underlying difference between the groups. For type II simulations, all group B words have slightly different targets (cf. Section 2.1 in the paper and also the Appendix).

The following code sets the parameters that determine the main characteristics of the data set.

``` r
# setting time dimension (these are the points along which the trajectories are generated)
xs_dense = seq(0,1,0.05)

# indices for extracting every second measurement; this is used for downsampling
xs_thin_ind = c(rep(c(T,F), (length(xs_dense)-1)/2), T) 

# population parameters: individual words come from this distribution
f2_start_mean = 1500
f2_end_1_mean = 2000
f2_end_2_mean = 1960

f2_start_sd.word = 40 # population level sd
f2_end_1_sd.word = 40 # population level sd
f2_end_2_sd.word = 40 # population level sd

# expected value & sd for transition point
x0_mean = 0.35
x0_sd.word = 0.020 # population level sd
# expected value & sd for steepness (higher -> more steep)
k_mean = 25
k_sd.word = 4 # population level sd

# how much variation within words?
f2_start_sd.traj = 30 # trajectory level sd
f2_end_1_sd.traj = 30 # trajectory level sd
f2_end_2_sd.traj = 30 # trajectory level sd
x0_sd.traj = 0.015 # trajectory level sd
k_sd.traj = 3 # trajectory level sd

# amount of random noise
noise_sd <- 5

# number of words & trajectories / word
n_words <- 50
n_trajectories_per_word <- 10
```

The code below assembles the dense version of the data set.

``` r
# creating matrix that will store the trajectories
ys_m <- matrix(0, nrow=length(xs_dense), ncol=n_words*n_trajectories_per_word)

# assembling words & individual trajectories
for (i in 1:n_words) {
  # sampling word-level targets
  f2_start.word <- rnorm(1, f2_start_mean, f2_start_sd.word)
  f2_end_1.word <- rnorm(1, f2_end_1_mean, f2_end_1_sd.word)
  if (type==1) {
    f2_end_2.word <- rnorm(1, f2_end_1_mean, f2_end_1_sd.word)
  } else {
    f2_end_2.word <- rnorm(1, f2_end_2_mean, f2_end_2_sd.word)
  }
  x0.word <- rnorm(1, x0_mean, x0_sd.word)
  k.word <- rnorm(1, k_mean, k_sd.word)
  
  # sampling trajectory-level targets,
  for (j in 1:n_trajectories_per_word) {
    f2_start <- rnorm(1, f2_start.word, f2_start_sd.traj)
    f2_end_1 <- rnorm(1, f2_end_1.word, f2_end_1_sd.traj)
    f2_end_2 <- rnorm(1, f2_end_2.word, f2_end_2_sd.traj)
    x0 <- rnorm(1, x0.word, x0_sd.traj)
    k <- rnorm(1, k.word, k_sd.traj)
    
    # assembling trajectories
    if (i <= (n_words / 2)) {
      ys_m[,(i-1)*n_trajectories_per_word + j] <- ((f2_end_1 - f2_start) / (1 + exp(-k*(xs_dense-x0)))) + f2_start + rnorm(length(xs_dense), 0, noise_sd)
    } else {
      ys_m[,(i-1)*n_trajectories_per_word + j] <- ((f2_end_2 - f2_start) / (1 + exp(-k*(xs_dense-x0)))) + f2_start + rnorm(length(xs_dense), 0, noise_sd)
    }
  }
}

# assembling data set (randomly assigned to categories)
dat_dense <- data.frame(traj=paste("traj_", rep(1:(n_words*n_trajectories_per_word), each=length(xs_dense)), sep=""),
                  word=paste("word_", rep(1:n_words, each=length(xs_dense)*n_trajectories_per_word), sep=""),
                  group=rep(c("A","B"), each=length(xs_dense)*(n_words*n_trajectories_per_word / 2)),
                  measurement.no=xs_dense, 
                  f2=c(ys_m),
                  stringsAsFactors = F
                 )

# setting up different types of grouping factors
dat_dense$group.factor <- as.factor(dat_dense$group)
dat_dense$group.ordered <- as.ordered(dat_dense$group)
contrasts(dat_dense$group.ordered) <- "contr.treatment"
dat_dense$group.bin <- as.numeric(dat_dense$group.factor) - 1

# traj/word ids must be factors  
dat_dense$traj <- as.factor(dat_dense$traj)
dat_dense$word <- as.factor(dat_dense$word)

# add dat$start for AR.start (for autoregressive error models)

dat_dense$start <- dat_dense$measurement.no == 0
```

Downsampling to thin version of the data set.

``` r
dat_thin <- dat_dense[rep(xs_thin_ind, n_words*n_trajectories_per_word),]
```

Here is what the data set looks like.

``` r
ggplot(dat_dense, aes(x=measurement.no, y=f2, group=traj, col=group)) +
  geom_line() +
  facet_wrap(~word)
```

![](set_2-f2_files/figure-markdown_github/unnamed-chunk-6-1.png)

Models
------

All the models (and sets of models) from Table 3 are shown below in the same order as in the table. Note that all models contain AR1 components to deal with dependencies within trajectories. For simplicity, the rho value used for these AR1 components is taken from a single model fitted without any random structures. This model is estimated below.

``` r
# dense
rho_mod_dense <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10), 
                    data = dat_dense, method = "fREML", discrete = T, nthreads = 1)

# thin
rho_mod_thin <- bam(f2 ~ group.ordered + 
                     s(measurement.no, bs = "tp", k = 10) + 
                     s(measurement.no, by = group.ordered, bs = "tp", k = 10), 
                   data = dat_thin, method = "fREML", discrete = T, nthreads = 1)

rho_dense <- start_value_rho(rho_mod_dense)
rho_thin <- start_value_rho(rho_mod_thin)
```

### 1. No components

``` r
# dense
nocomp_dense <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10), 
                    data = dat_dense, 
                    AR.start = dat_dense$start, rho = rho_dense, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(nocomp_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1814.890      2.066 878.336  < 2e-16 ***
    ## group.orderedB  -15.795      2.921  -5.406 6.57e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df       F p-value    
    ## s(measurement.no)                8.994  8.999 5044.26  <2e-16 ***
    ## s(measurement.no):group.orderedB 7.493  8.559   30.86  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.932   Deviance explained = 93.2%
    ## fREML =  45864  Scale est. = 1892.5    n = 10500

``` r
# thin
nocomp_thin <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10), 
                    data = dat_thin, 
                    AR.start = dat_thin$start, rho = rho_thin, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(nocomp_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1811.525      2.099 863.068  < 2e-16 ***
    ## group.orderedB  -15.645      2.967  -5.272  1.4e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df       F p-value    
    ## s(measurement.no)                8.990  8.998 3647.81  <2e-16 ***
    ## s(measurement.no):group.orderedB 6.609  7.932   23.95  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.934   Deviance explained = 93.5%
    ## fREML =  26767  Scale est. = 2277      n = 5500

### 2. Rand intcpt (= Random Intercept)

``` r
# dense
rand_intcpt_dense <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10) +
                      s(word, bs = "re"), 
                    data = dat_dense, 
                    AR.start = dat_dense$start, rho = rho_dense, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(word, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1794.179      5.760 311.470   <2e-16 ***
    ## group.orderedB  -13.778      8.146  -1.691   0.0908 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df        F p-value    
    ## s(measurement.no)                 8.994  8.999 5220.791  <2e-16 ***
    ## s(measurement.no):group.orderedB  7.534  8.582   31.909  <2e-16 ***
    ## s(word)                          42.548 48.000    7.805  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.952   Deviance explained = 95.3%
    ## fREML =  45732  Scale est. = 1827.2    n = 10500

``` r
# thin
rand_intcpt_thin <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10) +
                      s(word, bs = "re"), 
                    data = dat_thin, 
                    AR.start = dat_thin$start, rho = rho_thin, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(word, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1833.284      5.839 313.995   <2e-16 ***
    ## group.orderedB  -17.531      8.256  -2.123   0.0338 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df        F p-value    
    ## s(measurement.no)                 8.990  8.999 3890.885  <2e-16 ***
    ## s(measurement.no):group.orderedB  6.707  8.009   25.417  <2e-16 ***
    ## s(word)                          42.568 48.000    7.836  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.954   Deviance explained = 95.4%
    ## fREML =  26638  Scale est. = 2130.8    n = 5500

### 3. Rand intcpt + slope (= random intercept + slope)

``` r
# dense
rand_intcpt_slope_dense <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10) +
                      s(word, bs = "re") +
                      s(word, measurement.no, bs = "re"), 
                    data = dat_dense, 
                    AR.start = dat_dense$start, rho = rho_dense, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_slope_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(word, bs = "re") + 
    ##     s(word, measurement.no, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1794.175      8.824 203.321   <2e-16 ***
    ## group.orderedB  -13.771     12.479  -1.104     0.27    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df       F p-value    
    ## s(measurement.no)                 8.994  8.999 3075.20  <2e-16 ***
    ## s(measurement.no):group.orderedB  7.598  8.616   14.98  <2e-16 ***
    ## s(word)                          41.649 48.000  109.74  <2e-16 ***
    ## s(measurement.no,word)           44.450 48.000  126.33  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.964   Deviance explained = 96.4%
    ## fREML =  45503  Scale est. = 1725.7    n = 10500

``` r
# thin
rand_intcpt_slope_thin <- bam(f2 ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 10) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 10) +
                      s(word, bs = "re") +
                      s(word, measurement.no, bs = "re"), 
                    data = dat_thin, 
                    AR.start = dat_thin$start, rho = rho_thin, 
                    method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_slope_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(word, bs = "re") + 
    ##     s(word, measurement.no, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1833.294      8.882 206.412   <2e-16 ***
    ## group.orderedB  -17.549     12.560  -1.397    0.162    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df       F p-value    
    ## s(measurement.no)                 8.991  8.999 2187.56 < 2e-16 ***
    ## s(measurement.no):group.orderedB  6.838  8.109   10.83 2.3e-15 ***
    ## s(word)                          41.214 48.000  110.75 < 2e-16 ***
    ## s(measurement.no,word)           43.864 48.000  125.28 < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.966   Deviance explained = 96.6%
    ## fREML =  26455  Scale est. = 1944.6    n = 5500

### 4. Rand smooth, 3 bs

``` r
# dense
rand_smooth_3_dense <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 3), 
      data = dat_dense, 
      AR.start = dat_dense$start, rho = rho_dense, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_3_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 3)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1788.897      6.242  286.60   <2e-16 ***
    ## group.orderedB  -11.605      8.723   -1.33    0.183    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df        F  p-value    
    ## s(measurement.no)                  8.994   8.999 1692.081  < 2e-16 ***
    ## s(measurement.no):group.orderedB   7.411   8.405    6.486 1.18e-08 ***
    ## s(measurement.no,word)           131.478 148.000   11.880  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =   0.97   Deviance explained = 97.1%
    ## fREML =  45232  Scale est. = 1620.7    n = 10500

``` r
# thin
rand_smooth_3_thin <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 3), 
      data = dat_thin, 
      AR.start = dat_thin$start, rho = rho_thin, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_3_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 3)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1831.126      6.367 287.593   <2e-16 ***
    ## group.orderedB  -15.442      8.905  -1.734    0.083 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df        F  p-value    
    ## s(measurement.no)                  8.991   8.998 1318.856  < 2e-16 ***
    ## s(measurement.no):group.orderedB   6.731   7.969    5.142 2.18e-06 ***
    ## s(measurement.no,word)           128.663 148.000    9.844  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.971   Deviance explained = 97.2%
    ## fREML =  26288  Scale est. = 1799      n = 5500

### 5. Rand smooth, 5 bs

``` r
# dense
rand_smooth_5_dense <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 5), 
      data = dat_dense, 
      AR.start = dat_dense$start, rho = rho_dense, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_5_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 5)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1779.763      6.811 261.293   <2e-16 ***
    ## group.orderedB  -13.801      9.161  -1.506    0.132    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df       F p-value    
    ## s(measurement.no)                  8.991   8.995 931.543  <2e-16 ***
    ## s(measurement.no):group.orderedB   5.250   6.266   1.738   0.113    
    ## s(measurement.no,word)           228.275 248.000  16.273  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.975   Deviance explained = 97.5%
    ## fREML =  44516  Scale est. = 1367.2    n = 10500

``` r
# thin
rand_smooth_5_thin <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 5), 
      data = dat_thin, 
      AR.start = dat_thin$start, rho = rho_thin, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_5_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 5)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1829.857      6.755 270.881   <2e-16 ***
    ## group.orderedB  -16.791      9.183  -1.828   0.0675 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df       F p-value    
    ## s(measurement.no)                  8.991   8.996 842.918  <2e-16 ***
    ## s(measurement.no):group.orderedB   4.248   5.297   2.113  0.0623 .  
    ## s(measurement.no,word)           222.020 248.000  11.968  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.975   Deviance explained = 97.6%
    ## fREML =  25894  Scale est. = 1479.2    n = 5500

### 6. Rand smooth, 10 bs

``` r
# dense
rand_smooth_10_dense <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 10), 
      data = dat_dense, 
      AR.start = dat_dense$start, rho = rho_dense, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_10_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 10)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1782.057      6.818 261.383   <2e-16 ***
    ## group.orderedB  -11.809      8.997  -1.313    0.189    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf  Ref.df       F p-value    
    ## s(measurement.no)                  8.96   8.965 367.620  <2e-16 ***
    ## s(measurement.no):group.orderedB   2.36   2.586   1.918   0.111    
    ## s(measurement.no,word)           471.02 498.000  22.450  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.977   Deviance explained = 97.8%
    ## fREML =  42883  Scale est. = 916.98    n = 10500

``` r
# thin
rand_smooth_10_thin <- 
  bam(f2 ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 10) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
        s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = 10), 
      data = dat_thin, 
      AR.start = dat_thin$start, rho = rho_thin, 
      method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(rand_smooth_10_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = 10)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1823.512      6.954 262.240   <2e-16 ***
    ## group.orderedB  -15.705      9.528  -1.648   0.0993 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df      F p-value    
    ## s(measurement.no)                  8.952   8.956 426.20  <2e-16 ***
    ## s(measurement.no):group.orderedB   2.767   2.830   3.10  0.0223 *  
    ## s(measurement.no,word)           457.173 496.000  14.69  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.976   Deviance explained = 97.8%
    ## fREML =  25144  Scale est. = 978.59    n = 5500

### 7. Rand smooth, gam.check

The model fitted in this section uses random smooths where the number of basis functions (k) is determined using the gam.check() function: after fitting an initial model with a relatively low value of k, gam.check() is used to see whether more wiggliness is necessary (essentially, whether the smooths use up all the degrees of freedom afforded to them). If gam.check() suggests that more wiggliness is necessary, this procedure is repeated again using a model with a higher value of k.

Below is a convenience function for extracting the relevant p-value from the output of gam.check.

``` r
gam.check.p.value <- function (mod, which.line) { # which.line is a regexp
  str.out <- capture.output(gam.check(mod))
  relevant.line <- str.out[grep(which.line, str.out)]
  p.value <- as.numeric(gsub(".*?([0-9.e-]*)[ .*]*$", "\\1", relevant.line, perl=T))
  return(p.value)
}
```

Fitting the models. Dense first.

``` r
# what k's should be tried?
k_min = 4
k_max = 10
k_step = 3
```

``` r
# dense
for (k in seq(k_min,k_max,k_step)) {
  cat("fitting model with  k =", k, "\n")
  
  # fitting model
  
  rand_smooth_gam.check_dense <- 
    bam(f2 ~ group.ordered + 
          s(measurement.no, bs = "tp", k = 10) + 
          s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
          s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = k), 
        data = dat_dense, 
        AR.start = dat_dense$start, rho = rho_dense, 
        method = "fREML", discrete = T, nthreads = 1)
  
  # check whether more complexity is needed using gam.check
  
  if (gam.check.p.value(rand_smooth_gam.check_dense, "word") >= 0.05 | k == k_max) {
    print(summary(rand_smooth_gam.check_dense))
    break
  }
}
```

    ## fitting model with  k = 4

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = k)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1786.734      6.504 274.729   <2e-16 ***
    ## group.orderedB  -11.128      9.019  -1.234    0.217    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df        F  p-value    
    ## s(measurement.no)                  8.993   8.997 1686.680  < 2e-16 ***
    ## s(measurement.no):group.orderedB   7.070   7.985    5.669 3.36e-07 ***
    ## s(measurement.no,word)           179.225 198.000   14.378  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.973   Deviance explained = 97.4%
    ## fREML =  44870  Scale est. = 1488.4    n = 10500

And now the thin data set.

``` r
# thin
for (k in seq(k_min,k_max,k_step)) {
  cat("fitting model with  k =", k, "\n")
  
  # fitting model
  
  rand_smooth_gam.check_thin <- 
    bam(f2 ~ group.ordered + 
          s(measurement.no, bs = "tp", k = 10) + 
          s(measurement.no, by = group.ordered, bs = "tp", k = 10) + 
          s(measurement.no, word, bs = "fs", m = 1, xt = "tp", k = k), 
        data = dat_thin, 
        AR.start = dat_thin$start, rho = rho_thin, 
        method = "fREML", discrete = T, nthreads = 1)
  
  # check whether more complexity is needed using gam.check
  
  if (gam.check.p.value(rand_smooth_gam.check_thin, "word") >= 0.05 | k == k_max) {
    print(summary(rand_smooth_gam.check_thin))
    break
  }
}
```

    ## fitting model with  k = 4

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f2 ~ group.ordered + s(measurement.no, bs = "tp", k = 10) + s(measurement.no, 
    ##     by = group.ordered, bs = "tp", k = 10) + s(measurement.no, 
    ##     word, bs = "fs", m = 1, xt = "tp", k = k)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    1830.300      6.624 276.321   <2e-16 ***
    ## group.orderedB  -14.771      9.206  -1.605    0.109    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df        F  p-value    
    ## s(measurement.no)                  8.991   8.997 1300.541  < 2e-16 ***
    ## s(measurement.no):group.orderedB   6.451   7.616    4.413 3.94e-05 ***
    ## s(measurement.no,word)           175.091 198.000   11.050  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.973   Deviance explained = 97.4%
    ## fREML =  26086  Scale est. = 1627.6    n = 5500
