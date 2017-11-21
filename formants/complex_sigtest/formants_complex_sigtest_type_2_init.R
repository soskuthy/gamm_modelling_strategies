#-----------------------------
# code for creating data set 
#-----------------------------

# type of data: simulated f2 trajectories for /eI/ (logistic curve)
# 500 trajectories, 50 words
# parameterised variation in
#   - start of trajectory
#   - end of trajectory
#   - transition point
#   - transition steepness
# + a bit of random noise
#   - hierarchical structure (measurements within trajectories within words)
#   - across-word variation

# setting time dimension

xs_dense = seq(0,1,0.05)
xs_thin_ind = c(rep(c(T,F), (length(xs_dense)-1)/2), T)

# population parameters: individual words come from this dist
f2_start_mean = 1500
f2_end_1_mean = 2000
f2_end_2_mean = 1960

f2_start_sd.word = 40
f2_end_1_sd.word = 40
f2_end_2_sd.word = 40
# expected value & sd for transition point
x0_mean = 0.35
x0_sd.word = 0.020
# expected value & sd for steepness (higher -> more steep)
k_mean = 25
k_sd.word = 4

# how much variation within words?
f2_start_sd.traj = 30
f2_end_1_sd.traj = 30
f2_end_2_sd.traj = 30
x0_sd.traj = 0.015
k_sd.traj = 3

# amount of random noise

noise_sd <- 5

n_words <- 50
n_trajectories_per_word <- 10

# assembling trajectories

ys_m <- matrix(0, nrow=length(xs_dense), ncol=n_words*n_trajectories_per_word)
for (i in 1:n_words) {
  f2_start.word <- rnorm(1, f2_start_mean, f2_start_sd.word)
  f2_end_1.word <- rnorm(1, f2_end_1_mean, f2_end_1_sd.word)
  f2_end_2.word <- rnorm(1, f2_end_2_mean, f2_end_2_sd.word)
  x0.word <- rnorm(1, x0_mean, x0_sd.word)
  k.word <- rnorm(1, k_mean, k_sd.word)
  for (j in 1:n_trajectories_per_word) {
    f2_start <- rnorm(1, f2_start.word, f2_start_sd.traj)
    f2_end_1 <- rnorm(1, f2_end_1.word, f2_end_1_sd.traj)
    f2_end_2 <- rnorm(1, f2_end_2.word, f2_end_2_sd.traj)
    x0 <- rnorm(1, x0.word, x0_sd.traj)
    k <- rnorm(1, k.word, k_sd.traj)
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
  fixef_formula <- "f2 ~ "
  nested_formula <- "f2 ~ "
  
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
    final_formula <- gsub("~ group.ordered", '~ s(group.ordered, bs="re")', final_formula)
    method_str <- 'method="fREML", select=TRUE'
  } else if (method == "discrete+select") {
    final_formula <- gsub("~ group.ordered", '~ s(group.ordered, bs="re")', final_formula)
    method_str <- 'method="fREML", discrete=T, nthreads=1, select=TRUE'
  }
  
  # assembling bam command
  bam_str <- paste("bam(", final_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat_", dataset, ", ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}
