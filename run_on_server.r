setwd("~/documents/research/projects/dynamic-gam/gamm_modelling_strategies")

# settings for current simulation
source("pitch/complex_sigtest/pitch_complex_sigtest_type_2_config_file.R")
min.t = 10101
max.t = 10240
time = "00:20:00"
mem = "1G"

config.file.path.on.server = file.path(path.on.server, config.file)
shell.script.name = paste(overall_name, ".sh", sep="")
job.file.name = paste(overall_name, ".job", sep="")
tar.on.server = paste(file.path(path.on.server, basename(output.dir)), '.tar.gz', sep="")
tar.on.local = paste(file.path(path.on.local, basename(output.dir)), '.tar.gz', sep="")
output.dir.on.local = file.path(path.on.local, basename(output.dir))
# shell script to execute on remote server:


shell.script = paste("#!/bin/sh
module load R/3.1.2
qsub ", path.on.server, "/", job.file.name, sep="")

# save shell script (set permissions?)

cat(shell.script, file=file.path(path.on.local, shell.script.name))
system(paste("chmod a=rwx ", file.path(path.on.local, shell.script.name)))
# job file

job.details = paste("#$ -cwd -V
#$ -l h_rt=", time, "
#$ -o ", path.on.server, "/logs
#$ -e ", path.on.server, "/logs
#$ -l h_vmem=", mem, "
#$ -t ", min.t, "-", max.t, "
#$ -N ", sep="")
job.name = paste(overall_name, "\n", sep="")
job.script = "Rscript --no-save --no-restore"

job.r.config.file.path = paste(config.file.path.on.server, "${SGE_TASK_ID}\n")

job.file = paste(job.details, job.name, "\n", job.script, " ", r.script.path.on.server, " ", job.r.config.file.path, sep="")
cat(job.file, file=file.path(path.on.local, job.file.name))

# files to copy to remote server (if not already there)

files.to.copy = c(file.path(path.on.local, shell.script.name),
                  file.path(path.on.local, job.file.name),
                  file.path(path.on.local, config.file),
                  file.path(path.on.local, init.file))


library(RCurl)
# refresh gamm single iteration:
# system(paste("export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; scp ", "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/gamm_single_iteration.r", " ms1341@login.yarcc.york.ac.uk:", r.script.path.on.server, sep=""))
# refresh library
# system(paste("export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; scp ", "~/documents/research/projects/dynamic-gam/gamm_modelling_strategies/libraries.r", " ms1341@login.yarcc.york.ac.uk:", file.path(dirname(r.script.path.on.server), "libraries.r"), sep=""))
system(paste('export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; ssh -l ms1341 login.yarcc.york.ac.uk "mkdir -p ', output.dir, '"', sep=""), wait=T)
system(paste("export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; scp ", paste(files.to.copy, collapse=" "), " ms1341@login.yarcc.york.ac.uk:", path.on.server, "/", sep=""))
system(paste('export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; ssh ms1341@login.yarcc.york.ac.uk "source /etc/profile; source ~/.bashrc; ', 'cd ', path.on.server, "; ", "./", shell.script.name, '"', sep=""))

# retrieve files
system(paste('export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; ssh ms1341@login.yarcc.york.ac.uk "tar -PC ', output.dir, ' -czf ', tar.on.server, ' . "', sep=""))
system(paste("export SSH_ASKPASS=/usr/local/Cellar/ssh-askpass/1.2.1/bin/ssh-askpass ; scp -r ms1341@login.yarcc.york.ac.uk:", tar.on.server, " ", tar.on.local, sep=""))
system(paste('mkdir -p ', output.dir.on.local, '; tar -xzf ', tar.on.local, ' -C ', output.dir.on.local, '/', " ; rm ", tar.on.local, sep=""))
