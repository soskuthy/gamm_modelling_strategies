# name and paths
overall_name = "gamm_formants_complex_across_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/formants/complex_across"
path.on.server = "/scratch/ms1341/r/gamm-sim/formants/complex_across"
output.dir = "/scratch/ms1341/r/gamm-sim/formants/complex_across/output_type_1"
r.script.path.on.server = "/scratch/ms1341/r/gamm-sim/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "formants_complex_across_type_1_init.r"

# name of config file
config.file = "formants_complex_across_type_1_config_file.R"

# details of models to be run 
to_fit <- expand.grid(fixed_effects=c("diff_tp_10"),
                      random_effects=c("noranef+noranef","noranef+rintcpt","noranef+rslope",
                                       "noranef+rsmooth_tp_3","noranef+rsmooth_tp_5","noranef+rsmooth_tp_10"),
                      AR=c("AR_est"),
                      method="discrete",
                      mod_comp="nomodcomp",
                      dataset=c("dense","thin"),
                      visual="noVis")


