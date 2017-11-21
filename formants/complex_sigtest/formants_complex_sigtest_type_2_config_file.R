# name and paths
overall_name = "gamm_formants_complex_sigtest_type_2"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/formants/complex_sigtest"
path.on.server = "/scratch/ms1341/r/gamm-sim/formants/complex_sigtest"
output.dir = "/scratch/ms1341/r/gamm-sim/formants/complex_sigtest/output_type_2"
r.script.path.on.server = "/scratch/ms1341/r/gamm-sim/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "formants_complex_sigtest_type_2_init.r"

# name of config file
config.file = "formants_complex_sigtest_type_2_config_file.R"

# details of models to be run 
to_fit <- rbind(expand.grid(fixed_effects=c("diff_tp_10"),
                            random_effects=c("noranef+rsmooth_tp_10"),
                            AR=c("AR_est"),
                            method=c("discrete","ML"),
                            mod_comp="modcomp",
                            dataset=c("thin"),
                            visual="vis"),
                expand.grid(fixed_effects=c("bin_tp_10"),
                            random_effects=c("noranef+rsmooth_tp_10"),
                            AR=c("AR_est"),
                            method=c("discrete"),
                            mod_comp="nomodcomp",
                            dataset=c("thin"),
                            visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_10"),
                            random_effects=c("noranef+rsmooth_tp_10"),
                            AR=c("AR_est"),
                            method=c("discrete+select"),
                            mod_comp="modcomp",
                            dataset=c("thin"),
                            visual="noVis")
)