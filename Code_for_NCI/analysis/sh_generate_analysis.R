project <- "vk72"
ncpus <- 48
mem_gb <- 192
walltime <- "24:00:00"

scripts <- c(manifest = "build_manifest.R", estimates = "build_estimates.R")

make_job <- function(script) {
  job <- sprintf('#!/bin/bash
#PBS -P %s
#PBS -q normal
#PBS -l ncpus=%d
#PBS -l mem=%dGB
#PBS -l walltime=%s
#PBS -l wd
#PBS -l storage=scratch/%s+gdata/%s

module load R/4.5.0

Rscript %s
', project, ncpus, mem_gb, walltime, project, project, script)
  gsub("\r", "", job)
}

for (nm in names(scripts)) {
  con <- file(sprintf("job_%s.sh", nm), open = "wb")
  writeLines(make_job(scripts[[nm]]), con, sep = "\n")
  close(con)
}

cat("wrote", paste0("job_", names(scripts), ".sh"), "\n")
