# name and paths
overall_name = "gamm_pitch_complex_sigtest_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/pitch/complex_sigtest"
path.on.server = "/scratch/ms1341/r/gamm-sim/pitch/complex_sigtest"
output.dir = "/scratch/ms1341/r/gamm-sim/pitch/complex_sigtest/output_type_1"
r.script.path.on.server = "/scratch/ms1341/r/gamm-sim/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "pitch_complex_sigtest_type_1_init.r"

# name of config file
config.file = "pitch_complex_sigtest_type_1_config_file.R"

# details of models to be run
to_fit <- rbind(expand.grid(fixed_effects=c("diff_tp_15"),
                            random_effects=c("noranef+rsmooth_tp_12"),
                            AR=c("AR_est"),
                            method=c("discrete","ML"),
                            mod_comp="modcomp",
                            dataset=c("dense"),
                            visual="vis"),
                expand.grid(fixed_effects=c("bin_tp_15"),
                            random_effects=c("noranef+rsmooth_tp_12"),
                            AR=c("AR_est"),
                            method=c("discrete"),
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_15"),
                            random_effects=c("noranef+rsmooth_tp_12"),
                            AR=c("AR_est"),
                            method=c("discrete+select"),
                            mod_comp="modcomp",
                            dataset=c("dense"),
                            visual="noVis"))


