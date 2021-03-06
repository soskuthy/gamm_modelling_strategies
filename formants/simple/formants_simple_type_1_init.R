#-----------------------------
# code for loading / resampling data set 
#-----------------------------

# type of data: resampled f2 trajectories for /aI/ (rising contour, simple shape)
# each simulation starts with contours from a different speaker
# each speaker has at least 20 trajectories! (10 per group)

# path previously determined via config file

dat_full <- readRDS(file.path(dirname(config.file.curr), data.file))

# limit to speakers with >= 20 trajectories

dat_full <- subset(dat_full, n >= 30)

# we sample a single speaker from this data set

dat <- subset(dat_full, speaker == sample(unique(dat_full$speaker), 1))

# we only take what we need
# (< 50 trajectories)

dat <- subset(dat, 
              traj %in% sample(unique(dat$traj), 
                               size=min(length(unique(dat$traj)), 50)
              )
)

# we now add randomly assigned category labels

ids <- unique(dat$traj)
group.Bs <- sample(ids, round(length(ids)/2))
dat$group <- "A"
dat$group[dat$traj %in% group.Bs] <- "B"

# setting up different types of grouping factors
dat$group.factor <- as.factor(dat$group)
dat$group.ordered <- as.ordered(dat$group)
contrasts(dat$group.ordered) <- "contr.treatment"
dat$group.bin <- as.numeric(dat$group.factor) - 1

# ids ought to be factors  
dat$traj <- as.factor(dat$traj)

# dat$start has already been added at data prep stage (for AR.start, i.e. for autoregressive error models)



#-----------------------------
# function for assembling bam
#-----------------------------

assemble_bam <- function (fixefs, ranefs, AR, method, dataset) {
  
  # extract parameters from fixef and ranef string specification
  fixef_specs <- str_split(fixefs, "_")[[1]]
  ranef_specs <- str_split(ranefs, "_")[[1]]
  
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
  } else if (ranef_specs[1] == "gamcheck") {
    ranef_formula <- paste(ranef_formula, 
                           '+ s(measurement.no, traj, bs="fs", m=1, xt="',
                           ranef_specs[2],
                           '", k=',
                           as.character(seq(as.numeric(ranef_specs[3]),as.numeric(ranef_specs[4]),as.numeric(ranef_specs[5]))),
                           ')',
                           sep="")
  } else if (ranef_specs[1] == "noranef") {
  } else {
    stop("Invalid random effect string: has to start with one of 'gamcheck', rintcpt', 'rslope', 'rsmooth' or 'noranef'")
  }
  
  final_formula <- paste(fixef_formula, ranef_formula)
  final_nested_formula <- paste(nested_formula, ranef_formula)
  
  # setting AR parameters based on string from job file
  
  AR_str <- ""
  if (AR == "AR_est") {
    AR_str <- paste("AR.start=dat$start, rho=rho.est, ", sep="")
  } else if (AR != "noAR") {
    AR_str <- paste("AR.start=dat$start, rho=", str_split(AR, "_")[[1]][2], sep="")
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
  bam_str <- paste("bam(", final_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}


#------------------------------------------------
# Function for extracting p-value from gam.check
#------------------------------------------------

gam.check.p.value <- function (mod, which.line) {
  str.out <- capture.output(gam.check(mod))
  relevant.line <- str.out[grep(which.line, str.out)]
  p.value <- as.numeric(str_match(relevant.line, "([0-9.]*)[ *.]*$")[[2]])
  return(p.value)
}
