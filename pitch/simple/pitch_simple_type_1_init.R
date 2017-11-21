# type of data: simulated pitch trajectories (in Hz)
# 50 trajectories, randomly assigned to two groups
# parameterised variation in
#   - start of trajectory
#   - overall declination (linear)
#   - boundary 1
#   - boundary 2
#   - H* - horizontal
#   - H* - vertical
#   - H- - vertical

# + a bit of random noise

# setting time dimension

xs_dense = seq(0,1,0.025)
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

# pitch accent

H.star_vertical_mean = 15 # roughly the equivalent of ~ 1.5-2 semitones
H.star_vertical_sd = 6

H.star_horizontal_min = 0 # as a proportion of the duration of W1
H.star_horizontal_max = 0.25

H.star_bw = 0.12 # as a proportion of *overall* duration

# boundary tone

H.minus_vertical_mean = 8 # about half of H.star
H.minus_vertical_sd = 3

H.minus_bw = 0.08

# final boundary tone

L.percent_vertical = -20
L.percent_bw = 0.12

# that's all!

n_trajectories <- 50

# assembling trajectories

ys_m <- matrix(0, nrow=length(xs_dense), ncol=n_trajectories)
for (i in 1:n_trajectories) {
  start <- rnorm(1, start_mean, start_sd)
  slope <- rnorm(1, slope_mean, slope_sd)
  boundary.1 <- runif(1, boundary.1_min, boundary.1_max)
  boundary.2 <- runif(1, boundary.2_min, boundary.2_max)
  H.star_horizontal <- runif(1, H.star_horizontal_min, H.star_horizontal_max)
  H.star_vertical <- rnorm(1, H.star_vertical_mean, H.star_vertical_sd)
  H.minus_vertical <- rnorm(1, H.minus_vertical_mean, H.minus_vertical_sd)
  ys_m[,i] <- start + xs_dense*slope +  # declination
    exp(-((xs_dense - (boundary.1*H.star_horizontal))**2)/(2*H.star_bw**2)) * H.star_vertical + # 1st pitch accent
    exp(-((xs_dense - boundary.2)**2)/(2*H.minus_bw**2)) * H.minus_vertical + # boundary tone
    exp(-((xs_dense - 1)**2)/(2*L.percent_bw**2)) * L.percent_vertical # final boundary tone
}

# assembling data set (randomly assigned categories)
dat_dense <- data.frame(traj=paste("traj_", rep(1:n_trajectories, each=length(xs_dense)), sep=""), 
                  group=rep(c("A","B"), each=length(xs_dense)*(n_trajectories / 2)),
                  measurement.no=xs_dense, 
                  pitch=c(ys_m),
                  stringsAsFactors = F
)

# setting up different types of grouping factors
dat_dense$group.factor <- as.factor(dat_dense$group)
dat_dense$group.ordered <- as.ordered(dat_dense$group) 
contrasts(dat_dense$group.ordered) <- "contr.treatment"
dat_dense$group.bin <- as.numeric(dat_dense$group.factor) - 1

# ids ought to be factors  
dat_dense$traj <- as.factor(dat_dense$traj)

# add dat$start for AR.start (for autoregressive error models)

dat_dense$start <- dat_dense$measurement.no == 0

# thin data set:

dat_thin <- dat_dense[rep(xs_thin_ind, n_trajectories),]


#-----------------------------
# function for assembling bam
#-----------------------------

assemble_bam <- function (fixefs, ranefs, AR, method, dataset) {
  
  # extract parameters from fixef and ranef string specification
  fixef_specs <- str_split(fixefs, "_")[[1]]
  ranef_specs <- str_split(ranefs, "_")[[1]]
  
  # initialise fixed effect part of formula
  fixef_formula <- "pitch ~ "
  nested_formula <- "pitch ~ "
  
  if (fixef_specs[1] == "diff") {
    # name of grouping variable (part of output)
    grouping_var <- "group.ordered"
    
    # if a difference smooth is used, assemble relevant fixef formula
    # using bs and k specifications from fixefs
    fixef_formula <- paste(fixef_formula,
                           'group.ordered + s(measurement.no, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ') + s(measurement.no, by=group.ordered, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ')',
                           sep="")
    nested_formula <- paste(nested_formula,
                            's(measurement.no, bs="',
                            fixef_specs[2],
                            '", k=',
                            fixef_specs[3],
                            ')',
                            sep="")
  } else if (fixef_specs[1] == "bin") {
    # name of grouping variable (part of output)
    grouping_var <- "group.bin"
    
    # if a binary smooth is used, assemble relevant fixef formula
    # using bs and k specifications from fixefs (basically, no parametric term)
    fixef_formula <- paste(fixef_formula,
                           's(measurement.no, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ') + s(measurement.no, by=group.bin, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ')',
                           sep="")
    nested_formula <- paste(nested_formula,
                            's(measurement.no, bs="',
                            fixef_specs[2],
                            '", k=',
                            fixef_specs[3],
                            ')',
                            sep="")
  } else {
    stop("Invalid 'by' specification: has to be either 'diff' or 'bin'.")
  }
  
  # initialise random effect formula
  ranef_formula <- ""
  
  # code for generating formula for random intercept / random intercept + slope
  # and random smooth (with user-specified k and basis type)
  if (ranef_specs[1] == "rintcpt") {
    ranef_formula <- paste(ranef_formula, '+ s(traj, bs="re")')
  } else if (ranef_specs[1] == "rslope") {
    ranef_formula <- paste(ranef_formula, '+ s(traj, bs="re") + s(traj, measurement.no, bs="re")')
  } else if (ranef_specs[1] == "rsmooth") {
    ranef_formula <- paste(ranef_formula, 
                           '+ s(measurement.no, traj, bs="fs", m=1, xt="',
                           ranef_specs[2],
                           '", k=',
                           ranef_specs[3],
                           ')',
                           sep="")
  } else if (ranef_specs[1] == "noranef") {
  } else {
    stop("Invalid random effect string: has to start with one of 'rintcpt', 'rslope', 'rsmooth' or 'noranef'")
  }
  
  final_formula <- paste(fixef_formula, ranef_formula)
  final_nested_formula <- paste(nested_formula, ranef_formula)
  
  # setting AR parameters based on string from job file
  
  AR_str <- ""
  if (AR == "AR_est") {
    AR_str <- paste("AR.start=dat_", dataset, "$start, rho=rho.est, ", sep="")
  } else if (AR != "noAR") {
    AR_str <- paste("AR.start=dat_", dataset, "$start, rho=", str_split(AR, "_")[[1]][2], sep="")
  }
  
  # setting method parameter(s) based on string from job file
  
  method_str <- ""
  if (method == "ML") {
    method_str <- 'method="ML"'
  } else if (method == "REML") {
    method_str <- 'method="REML"'
  } else if (method == "fREML") {
    method_str <- 'method="fREML"'
  } else if (method == "discrete") {
    method_str <- 'method="fREML", discrete=T, nthreads=1'
  }
  
  # assembling bam command
  bam_str <- paste("bam(", final_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}
