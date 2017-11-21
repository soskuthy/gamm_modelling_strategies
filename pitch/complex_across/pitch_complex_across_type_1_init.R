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

xs_dense = seq(0,1,0.025)
xs_thin_ind = c(rep(c(T,F), (length(xs_dense)-1)/2), T)

# expected values & sd for starting and end points -- population level
start_mean = 170
slope_mean = -30
start_sd.word = 8
slope_sd.word = 6

# trajectory-level variation
start_sd.traj = 6
slope_sd.traj = 5

# boundaries between N1, N2 and N3 -- population level

boundary.1_min.word = (1/3) - 0.08
boundary.1_max.word = (1/3) + 0.08 # orig sd: 0.058 ~ 0.035 + 0.046

boundary.2_min.word = (2/3) - 0.08
boundary.2_max.word = (2/3) + 0.08

# trajectory-level variation

boundary.1_sd.traj = 0.035
boundary.2_sd.traj = 0.035
truncate_boundary_distr_at = 0.1 # trajectories vary around word-level boundary targets as ~ truncnorm with a = -0.1, b = 0.1

# pitch accent

H.star_vertical_mean = 15 # roughly the equivalent of ~ 1.5-2 semitones
H.star_vertical_sd.word = 5
H.star_vertical_sd.traj = 3

H.star_horizontal_min = 0 # as a proportion of the duration of W1
H.star_horizontal_max = 0.25

H.star_bw = 0.12 # as a proportion of *overall* duration

# boundary tone

H.minus_vertical_mean = 8 # about half of H.star
H.minus_vertical_sd.word = 2.4
H.minus_vertical_sd.traj = 1.4

H.minus_bw = 0.08

# final boundary tone

L.percent_vertical = -20
L.percent_bw = 0.12

# that's all!

n_trajectories_per_word <- 10
n_words <- 50

# assembling trajectories

ys_m <- matrix(0, nrow=length(xs_dense), ncol=n_words*n_trajectories_per_word)
for (i in 1:n_words) {
  start.word <- rnorm(1, start_mean, start_sd.word)
  slope.word <- rnorm(1, slope_mean, slope_sd.word)
  boundary.1.word <- runif(1, boundary.1_min.word, boundary.1_max.word)
  boundary.2.word <- runif(1, boundary.2_min.word, boundary.2_max.word)
  
  H.star_vertical.word <- rnorm(1, H.star_vertical_mean, H.star_vertical_sd.word)
  H.minus_vertical.word <- rnorm(1, H.minus_vertical_mean, H.minus_vertical_sd.word)
  
  for (j in 1:n_trajectories_per_word) {
    start <- rnorm(1, start.word, start_sd.traj)
    slope <- rnorm(1, slope.word, slope_sd.traj)
    boundary.1 <- rtruncnorm(1, 
                             a=boundary.1.word-truncate_boundary_distr_at, 
                             b=boundary.1.word+truncate_boundary_distr_at,
                             mean=boundary.1.word,
                             sd=boundary.1_sd.traj)
    boundary.2 <- rtruncnorm(1, 
                             a=boundary.2.word-truncate_boundary_distr_at, 
                             b=boundary.2.word+truncate_boundary_distr_at,
                             mean=boundary.2.word,
                             sd=boundary.2_sd.traj)
    H.star_vertical <- rnorm(1, H.star_vertical.word, H.star_vertical_sd.traj)
    H.minus_vertical <- rnorm(1, H.minus_vertical.word, H.minus_vertical_sd.traj)
    # this one varies randomly within trajs:  
    H.star_horizontal <- runif(1, H.star_horizontal_min, H.star_horizontal_max)

    ys_m[,(i-1)*n_trajectories_per_word + j] <- start + xs_dense*slope +  # declination
      exp(-((xs_dense - (boundary.1*H.star_horizontal))**2)/(2*H.star_bw**2)) * H.star_vertical + # 1st pitch accent
      exp(-((xs_dense - boundary.2)**2)/(2*H.minus_bw**2)) * H.minus_vertical + # boundary tone
      exp(-((xs_dense - 1)**2)/(2*L.percent_bw**2)) * L.percent_vertical # final boundary tone
  }
}

dat_dense <- data.frame(traj=paste("traj_", rep(1:(n_words*n_trajectories_per_word), each=length(xs_dense)), sep=""),
                        word=paste("word_", rep(1:n_words, each=length(xs_dense)*n_trajectories_per_word), sep=""),
                        group=rep(c("A","B"), each=length(xs_dense)*(n_words*n_trajectories_per_word / 2)),
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
dat_dense$word <- as.factor(dat_dense$word)

# add dat$start for AR.start (for autoregressive error models)

dat_dense$start <- dat_dense$measurement.no == 0

# thin data set:

dat_thin <- dat_dense[rep(xs_thin_ind, n_words*n_trajectories_per_word),]


#-----------------------------
# function for assembling bam
#-----------------------------

### utility function for assembling random structure
# INPUT: rintcpt OR rslope OR rsmooth_(tp/cr/..)_(k), e.g. rsmooth_cr_3
#        name of grouping variable
assemble_rstruct <- function (r, grouping) {
  out <- ''
  if (r[1] == "rintcpt") {
    out <- paste(out, ' + s(',  grouping, ', bs="re")', sep="")
  } else if (r[1] == "rslope") {
    out <- paste(out, ' + s(',  grouping, ', bs="re") + s(',  grouping, ', measurement.no, bs="re")', sep="")
  } else if (r[1] == "rsmooth") {
    out <- paste(out, 
                           ' + s(measurement.no, ',  grouping, ', bs="fs", m=1, xt="',
                           r[2],
                           '", k=',
                           r[3],
                           ')',
                           sep="")
  } else if (r[1] == "noranef") {
  } else {
    stop("Invalid random effect string: has to start with one of 'rintcpt', 'rslope', 'rsmooth' or 'noranef'")
  }
  return(out)
}

assemble_bam <- function (fixefs, ranefs, AR, method, dataset) {
  
  # extract parameters from fixef and ranef string specification
  fixef_specs <- str_split(fixefs, "_")[[1]]
  ranef_specs <- str_split(unlist(str_split(ranefs, "[+]")[[1]]), "_")
  
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
  
  # code for generating formula for random intercept / random intercept + slope
  # and random smooth (with user-specified k and basis type)
  
  ranef_formula <- ''
  groupers <- c("traj", "word")
  for (j in 1:length(ranef_specs)) {
    ranef_formula <- paste(ranef_formula, assemble_rstruct(ranef_specs[[j]], groupers[j]), sep="")
  }
  
  final_formula <- paste(fixef_formula, ranef_formula, sep="")
  final_nested_formula <- paste(nested_formula, ranef_formula, sep="")
  
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
  } else if (method == "fREML+select") {
    method_str <- 'method="fREML", select=TRUE'
  } else if (method == "discrete+select") {
    method_str <- 'method="fREML", discrete=T, nthreads=1, select=TRUE'
  }
  
  # assembling bam command
  bam_str <- paste("bam(", final_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}
