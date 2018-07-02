Introduction
------------

This document is a supplement to "Evaluating generalised additive mixed modelling strategies for dynamic speech analysis," relating specifically to the contents of Table 2 in Section 3.1. It presents code that illustrates (i) how the simulated data were generated and (ii) the models whose performance is summarised Table 1.

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

The code in this section can be used to create data for either type I or type II simulations. Set the value of *type* to 1 for type I simulations and to 2 for type II simulations. To replicate the results in the table in appendix C, which summarises the same types of trajectories, but with a small amount of added noise, set the value of *noise* to TRUE.

``` r
type = 1
noise = FALSE
```

The data for this set of simulations consist of simulated pitch trajectories loosely modelled after triconstituent compounds in English. 50 trajectories are generated. For type I simulations, these are randomly assigned to two groups (A and B). For type II simulations, all group B trajectories are slightly modified (cf. Section 2.1 in the paper and also the Appendix).

The following code sets the parameters that determine the main characteristics of the data set.

``` r
# setting time dimension (these are the points along which the trajectories are generated)
xs_dense = seq(0,1,0.025)

# indices for extracting every second measurement; this is used for downsampling
xs_thin_ind = c(rep(c(T,F), (length(xs_dense)-1)/2), T)

# expected values & sd for starting and end points
start_mean = 170
start_sd = 10
slope_mean = -30
slope_sd = 8

# boundaries between N1, N2 and N3

boundary.1_min = (1/3) - 0.1
boundary.1_max = (1/3) + 0.1

boundary.2_min = (2/3) - 0.1
boundary.2_max = (2/3) + 0.1

# pitch accent 1 - vertical target
H.star.1_vertical_mean = 15 # roughly the equivalent of ~ 1.5-2 semitones
H.star.1_vertical_sd = 6

# pitch accent 1 - horizontal target
H.star.1_horizontal_min = 0 # as a proportion of the duration of W1
H.star.1_horizontal_max = 0.25

# pitch accent 1 - width
H.star.1_bw = 0.12 # as a proportion of *overall* duration

# pitch accent 2 - vertical target
H.star.2_vertical_mean = 4.5
H.star.2_vertical_sd = 6

# pitch accent 2 - horizontal target
H.star.2_horizontal_min = 0 # as a proportion of the duration of W2
H.star.2_horizontal_max = 0.25

# pitch accent 2 - width
H.star.2_bw = 0.12 # as a proportion of *overall* duration
  
# boundary tone - vertical target
H.minus_vertical_mean = 8 # about half of H.star
H.minus_vertical_sd = 3

# boundary tone - width
H.minus_bw = 0.08

# final boundary tone - vertical target & width
L.percent_vertical = -20
L.percent_bw = 0.12

# noise
noise_sd = 1

# number of trajectories to generate
n_trajectories <- 50
```

The code below assembles the dense version of the data set.

``` r
# creating matrix that will store the trajectories
ys_m <- matrix(0, nrow=length(xs_dense), ncol=n_trajectories)

# assembling individual trajectories in pairs (one for group A and one for group B)
for (i in 1:n_trajectories) {
  # sampling start values, slopes, boundary positions and pitch accents / boundary tones 
  # for trajectories A and B
  start <- rnorm(1, start_mean, start_sd)
  slope <- rnorm(1, slope_mean, slope_sd)
  boundary.1 <- runif(1, boundary.1_min, boundary.1_max)
  boundary.2 <- runif(1, boundary.2_min, boundary.2_max)
  H.star.1_horizontal <- runif(1, H.star.1_horizontal_min, H.star.1_horizontal_max)
  H.star.1_vertical <- rnorm(1, H.star.1_vertical_mean, H.star.1_vertical_sd)
  H.star.2_horizontal <- runif(1, H.star.2_horizontal_min, H.star.2_horizontal_max)
  H.star.2_vertical <- rnorm(1, H.star.2_vertical_mean, H.star.2_vertical_sd)
  H.minus_vertical <- rnorm(1, H.minus_vertical_mean, H.minus_vertical_sd)
  
  # assembling the trajectories
  ys_m[,i] <- start + xs_dense*slope +  # declination
    exp(-((xs_dense - (boundary.1*H.star.1_horizontal))**2)/(2*H.star.1_bw**2)) * H.star.1_vertical + # 1st pitch accent
    exp(-((xs_dense - boundary.2)**2)/(2*H.minus_bw**2)) * H.minus_vertical + # boundary tone
    exp(-((xs_dense - 1)**2)/(2*L.percent_bw**2)) * L.percent_vertical # final boundary tone
  
  # adding additional pitch accent for group B if this is a type II simulation
  if (i > (n_trajectories/2) & type==2) {
    ys_m[,i] <- ys_m[,i] + 
      exp(-((xs_dense - (boundary.1+(boundary.2-boundary.1)*H.star.2_horizontal))**2)/(2*H.star.2_bw**2)) * H.star.2_vertical # 2nd pitch accent
  }
  
  # adding noise
  if (noise) {
    ys_m[,i] <- ys_m[,i] + rnorm(length(xs_dense), 0, noise_sd)
  }
}

# assembling data set
dat_dense <- data.frame(traj=paste("traj_", rep(1:n_trajectories, each=length(xs_dense)), sep=""), # traj ids
                  group=rep(c("A","B"), each=length(xs_dense)*(n_trajectories / 2)),               # group ids
                  measurement.no=xs_dense,                          # measurement number (i.e. time dimension)
                  pitch=c(ys_m),                                                             # pitch  values
                  stringsAsFactors = F
)

# setting up different types of grouping factors
dat_dense$group.factor <- as.factor(dat_dense$group)
dat_dense$group.ordered <- as.ordered(dat_dense$group) 
contrasts(dat_dense$group.ordered) <- "contr.treatment"
dat_dense$group.bin <- as.numeric(dat_dense$group.factor) - 1

# trajectory ids must be factors  
dat_dense$traj <- as.factor(dat_dense$traj)

# add dat$start for AR.start (for autoregressive error models)
dat_dense$start <- dat_dense$measurement.no == 0
```

Downsampling to thin version of the data set.

``` r
dat_thin <- dat_dense[rep(xs_thin_ind, n_trajectories),]
```

Here is what the data set looks like.

``` r
ggplot(dat_dense, aes(x=measurement.no, y=pitch, group=traj, col=group)) +
  geom_line() +
  facet_grid(~group)
```

![](set_1-pitch_files/figure-markdown_github/unnamed-chunk-6-1.png)

Models
------

All the models (and sets of models) from Table 1 are shown below in the same order as in the table.

### 1. No components

``` r
# dense
nocomp_dense <- bam(pitch ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 15) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 15),
                    data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
summary(nocomp_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    158.4112     0.3855 410.924   <2e-16 ***
    ## group.orderedB  -1.2999     0.5452  -2.384   0.0172 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df       F p-value    
    ## s(measurement.no)                8.02  9.748 229.180  <2e-16 ***
    ## s(measurement.no):group.orderedB 1.00  1.001   0.205   0.651    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.666   Deviance explained = 66.8%
    ## fREML = 8069.1  Scale est. = 152.33    n = 2050

``` r
# thin
nocomp_thin <- bam(pitch ~ group.ordered + 
                     s(measurement.no, bs = "tp", k = 15) + 
                     s(measurement.no, by = group.ordered, bs = "tp", k = 15), 
                   data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
summary(nocomp_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    158.3004     0.5413 292.446   <2e-16 ***
    ## group.orderedB  -1.3096     0.7655  -1.711   0.0874 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df       F p-value    
    ## s(measurement.no)                6.98  8.570 139.271  <2e-16 ***
    ## s(measurement.no):group.orderedB 1.00  1.001   0.096   0.757    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.677   Deviance explained =   68%
    ## fREML =   4139  Scale est. = 153.83    n = 1050

### 2. Rand intcpt (= Random Intercept)

``` r
# dense
rand_intcpt_dense <- bam(pitch ~ group.ordered + 
                           s(measurement.no, bs = "tp", k = 15) + 
                           s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
                           s(traj, bs = "re"), 
                         data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(traj, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.614      2.389  66.384   <2e-16 ***
    ## group.orderedB   -1.299      3.379  -0.385    0.701    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df        F p-value    
    ## s(measurement.no)                12.01  13.35 1726.798  <2e-16 ***
    ## s(measurement.no):group.orderedB  1.00   1.00    2.113   0.146    
    ## s(traj)                          47.88  49.00  385.841  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.968   Deviance explained = 96.9%
    ## fREML = 5839.7  Scale est. = 14.818    n = 2050

``` r
# thin
rand_intcpt_thin <- bam(pitch ~ group.ordered + 
                           s(measurement.no, bs = "tp", k = 15) + 
                           s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
                           s(traj, bs = "re"), 
                         data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
summary(rand_intcpt_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(traj, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.130      2.393  66.094   <2e-16 ***
    ## group.orderedB   -1.317      3.383  -0.389    0.697    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df       F p-value    
    ## s(measurement.no)                11.09  12.71 926.000  <2e-16 ***
    ## s(measurement.no):group.orderedB  1.00   1.00   0.942   0.332    
    ## s(traj)                          47.75  48.00 190.499  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.967   Deviance explained = 96.9%
    ## fREML = 3080.6  Scale est. = 15.686    n = 1050

### 3. Rand intcpt + slope (= random intercept + slope)

``` r
# dense
rand_intcpt_slope_dense <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(traj, bs = "re") + 
        s(traj, measurement.no, bs = "re"), 
      data = dat_dense, method = "fREML", discrete = T, nthreads = 1) 
summary(rand_intcpt_slope_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(traj, bs = "re") + s(traj, measurement.no, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.613      2.429  65.305   <2e-16 ***
    ## group.orderedB   -1.295      3.435  -0.377    0.706    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df         F p-value    
    ## s(measurement.no)                13.131 13.846 4.546e+02  <2e-16 ***
    ## s(measurement.no):group.orderedB  4.715  5.854 2.123e+00  0.0602 .  
    ## s(traj)                          47.822 48.000 2.579e+05  <2e-16 ***
    ## s(measurement.no,traj)           47.462 48.000 2.552e+05  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =   0.99   Deviance explained =   99%
    ## fREML = 4813.1  Scale est. = 4.7298    n = 2050

``` r
# thin
rand_intcpt_slope_thin <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(traj, bs = "re") + 
        s(traj, measurement.no, bs = "re"), 
      data = dat_thin, method = "fREML", discrete = T, nthreads = 1) 
summary(rand_intcpt_slope_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(traj, bs = "re") + s(traj, measurement.no, bs = "re")
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.122      2.415  65.469   <2e-16 ***
    ## group.orderedB   -1.315      3.415  -0.385      0.7    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df         F p-value    
    ## s(measurement.no)                12.604 13.664   392.318  <2e-16 ***
    ## s(measurement.no):group.orderedB  1.458  1.783     0.323   0.745    
    ## s(traj)                          47.643 48.000 57480.101  <2e-16 ***
    ## s(measurement.no,traj)           46.915 48.000 56279.503  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.989   Deviance explained =   99%
    ## fREML = 2609.9  Scale est. = 5.0721    n = 1050

### 4. Rand smooth, 5 bs

``` r
# dense
rand_smooth_3_dense <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 5), 
      data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 5)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.867      2.584  61.469   <2e-16 ***
    ## group.orderedB   -1.841      3.642  -0.506    0.613    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df       F p-value    
    ## s(measurement.no)                 13.685  13.94 161.369  <2e-16 ***
    ## s(measurement.no):group.orderedB   5.225   6.54   2.181  0.0663 .  
    ## s(measurement.no,traj)           237.517 249.00 965.746  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.997   Deviance explained = 99.8%
    ## fREML = 3825.6  Scale est. = 1.2793    n = 2050

``` r
# thin
rand_smooth_3_thin <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 5), 
      data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 5)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.374      2.571  61.605   <2e-16 ***
    ## group.orderedB   -1.838      3.617  -0.508    0.611    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df       F p-value    
    ## s(measurement.no)                 13.442  13.92 159.408  <2e-16 ***
    ## s(measurement.no):group.orderedB   2.503   3.05   0.155   0.937    
    ## s(measurement.no,traj)           232.203 249.00 405.288  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.997   Deviance explained = 99.7%
    ## fREML = 2273.2  Scale est. = 1.5659    n = 1050

### 5. Rand smooth, 8 bs

``` r
# dense
rand_smooth_5_dense <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 8), 
      data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 8)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.924      2.592  61.320   <2e-16 ***
    ## group.orderedB   -2.141      3.641  -0.588    0.557    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf  Ref.df        F p-value    
    ## s(measurement.no)                 13.85  13.893  180.019  <2e-16 ***
    ## s(measurement.no):group.orderedB   4.37   5.524    0.267   0.962    
    ## s(measurement.no,traj)           387.06 398.000 4794.094  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =      1   Deviance explained =  100%
    ## fREML = 2300.6  Scale est. = 0.16241   n = 2050

``` r
# thin
rand_smooth_5_thin <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 8), 
      data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 8)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.397      2.578  61.441   <2e-16 ***
    ## group.orderedB   -1.895      3.612  -0.525      0.6    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df        F p-value    
    ## s(measurement.no)                 13.822  13.906  184.289  <2e-16 ***
    ## s(measurement.no):group.orderedB   1.435   1.631    0.085   0.913    
    ## s(measurement.no,traj)           383.388 398.000 1739.768  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =      1   Deviance explained =  100%
    ## fREML = 1765.3  Scale est. = 0.23021   n = 1050

### 6. Rand smooth, 12 bs

``` r
# dense
rand_smooth_10_dense <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 12), 
      data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 12)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.397      2.557  61.939   <2e-16 ***
    ## group.orderedB   -1.910      3.588  -0.532    0.595    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df         F p-value    
    ## s(measurement.no)                 13.849  13.856   271.137  <2e-16 ***
    ## s(measurement.no):group.orderedB   2.009   2.482     2.044   0.151    
    ## s(measurement.no,traj)           587.893 599.000 61226.262  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =      1   Deviance explained =  100%
    ## fREML = 461.03  Scale est. = 0.0084581  n = 2050

``` r
# thin
rand_smooth_10_thin <- 
  bam(pitch ~ group.ordered + 
        s(measurement.no, bs = "tp", k = 15) + 
        s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
        s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 12), 
      data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = 12)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.096      2.554  61.899   <2e-16 ***
    ## group.orderedB   -1.726      3.581  -0.482     0.63    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                     edf Ref.df        F p-value    
    ## s(measurement.no)                 13.83  13.84   260.89  <2e-16 ***
    ## s(measurement.no):group.orderedB   1.00   1.00     0.72   0.397    
    ## s(measurement.no,traj)           586.12 599.00 25170.59  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =      1   Deviance explained =  100%
    ## fREML = 1264.9  Scale est. = 0.010589  n = 1050

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
k_max = 12
k_step = 4
```

``` r
# dense
for (k in seq(k_min,k_max,k_step)) {
  cat("fitting model with  k =", k, "\n")
  
  # fitting model
  
  rand_smooth_gam.check_dense <- 
    bam(pitch ~ group.ordered + 
          s(measurement.no, bs = "tp", k = 15) + 
          s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
          s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = k), 
        data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
  
  # check whether more complexity is needed using gam.check
  
  if (gam.check.p.value(rand_smooth_gam.check_dense, "traj") >= 0.05 | k == k_max) {
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = k)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.600      2.465  64.338   <2e-16 ***
    ## group.orderedB   -1.272      3.481  -0.365    0.715    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df       F  p-value    
    ## s(measurement.no)                 13.457  13.906 204.867  < 2e-16 ***
    ## s(measurement.no):group.orderedB   5.939   7.349   3.713 0.000378 ***
    ## s(measurement.no,traj)           187.625 199.000 616.663  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.995   Deviance explained = 99.5%
    ## fREML = 4341.7  Scale est. = 2.4868    n = 2050

And now the thin data set.

``` r
# thin
for (k in seq(k_min,k_max,k_step)) {
  cat("fitting model with  k =", k, "\n")
  
  # fitting model
  
  rand_smooth_gam.check_thin <- 
    bam(pitch ~ group.ordered + 
          s(measurement.no, bs = "tp", k = 15) + 
          s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
          s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = k), 
        data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
  
  # check whether more complexity is needed using gam.check
  
  if (gam.check.p.value(rand_smooth_gam.check_thin, "traj") >= 0.05 | k == k_max) {
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
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15) + 
    ##     s(measurement.no, traj, bs = "fs", m = 1, xt = "tp", k = k)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     158.075      2.465  64.140   <2e-16 ***
    ## group.orderedB   -1.215      3.479  -0.349    0.727    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df       F p-value    
    ## s(measurement.no)                 13.069  13.821 184.430  <2e-16 ***
    ## s(measurement.no):group.orderedB   4.027   5.007   1.512    0.19    
    ## s(measurement.no,traj)           182.022 200.000 274.659  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.994   Deviance explained = 99.5%
    ## fREML = 2454.6  Scale est. = 2.8524    n = 1050

### 8. AR1

First fitting models without AR component in order to estimate rho. This is equivalent to the nocomp model above.

``` r
# dense
nocomp_dense <- bam(pitch ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 15) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 15), 
                    data = dat_dense, method = "fREML", discrete = T, nthreads = 1)
# thin
nocomp_thin <- bam(pitch ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 15) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 15), 
                    data = dat_thin, method = "fREML", discrete = T, nthreads = 1)
```

Extracting rho.

``` r
rho_dense <- start_value_rho(nocomp_dense)
cat("rho =", rho_dense, "for the dense data set\n")
```

    ## rho = 0.9669434 for the dense data set

``` r
rho_thin <- start_value_rho(nocomp_thin)
cat("rho =", rho_thin, "for the thin data set\n")
```

    ## rho = 0.9323774 for the thin data set

Fitting models with AR1

``` r
# dense
AR1_dense <- bam(pitch ~ group.ordered + 
                   s(measurement.no, bs = "tp", k = 15) + 
                   s(measurement.no, by = group.ordered, bs = "tp", k = 15), 
                 data = dat_dense, 
                 AR.start = dat_dense$start, rho = rho_dense, 
                 method = "fREML", discrete = T, nthreads = 1)
summary(AR1_dense)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    158.5390     0.5993 264.557   <2e-16 ***
    ## group.orderedB  -1.5552     0.8251  -1.885   0.0596 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df      F p-value    
    ## s(measurement.no)                13.73  13.99 564.06  <2e-16 ***
    ## s(measurement.no):group.orderedB  1.00   1.00   0.03   0.863    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.666   Deviance explained = 66.8%
    ## fREML = 2925.3  Scale est. = 14.23     n = 2050

``` r
# thin
AR1_thin <- bam(pitch ~ group.ordered + 
                   s(measurement.no, bs = "tp", k = 15) + 
                   s(measurement.no, by = group.ordered, bs = "tp", k = 15), 
                 data = dat_thin, 
                 AR.start = dat_thin$start, rho = rho_thin, 
                 method = "fREML", discrete = T, nthreads = 1)
summary(AR1_thin)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## pitch ~ group.ordered + s(measurement.no, bs = "tp", k = 15) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 15)
    ## 
    ## Parametric coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    158.4194     0.8172 193.857   <2e-16 ***
    ## group.orderedB  -1.5505     1.1272  -1.376    0.169    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                    edf Ref.df       F p-value    
    ## s(measurement.no)                13.41 13.952 287.597  <2e-16 ***
    ## s(measurement.no):group.orderedB  1.00  1.001   0.016     0.9    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.676   Deviance explained =   68%
    ## fREML = 2224.2  Scale est. = 26.997    n = 1050
