# name and paths
overall_name = "gamm_formants_simple_type_2"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/formants/simple"
path.on.server = "/scratch/ms1341/r/gamm-sim/formants/simple"
output.dir = "/scratch/ms1341/r/gamm-sim/formants/simple/output_type_2"
r.script.path.on.server = "/scratch/ms1341/r/gamm-sim/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "formants_simple_type_2_init.r"

# name of config file
config.file = "formants_simple_type_2_config_file.R"

# details of models to be run 
to_fit <- rbind(
  expand.grid(fixed_effects=c("diff_tp_10"),
              random_effects=c("noranef","rintcpt","rslope",
                               "rsmooth_tp_3","rsmooth_tp_5","rsmooth_tp_10"),
              AR=c("noAR"),
              method="discrete",
              mod_comp="nomodcomp",
              dataset=c("dense","thin"),
              visual="noVis"),
  expand.grid(fixed_effects=c("diff_tp_10"),
              random_effects="noranef",
              AR="AR_est",
              method="discrete",
              mod_comp="nomodcomp",
              dataset=c("dense","thin"),
              visual="noVis")
)
