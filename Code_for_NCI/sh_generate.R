project <- "vk72"
ncpus <- 48
mem_gb <- 192
walltime <- "08:00:00"

JQ_grid <- list(c(50, 100), c(100, 200), c(200, 400))
seeds <- 1:100
methods <- c("wb", "corrected", "fhem", "glmmTMB")
K <- 2

cmds <- character(0)
for (jq in JQ_grid) {
  for (s in seeds) {
    for (m in methods) {
      cmds <- c(cmds, sprintf("Rscript run_one_cell.R %d %d %d %s %d", jq[1], jq[2], s, m, K))
    }
  }
}

con <- file("cmds.txt", open = "wb")
writeLines(cmds, con, sep = "\n")
close(con)

job_script <- sprintf('#!/bin/bash
#PBS -P %s
#PBS -q normal
#PBS -l ncpus=%d
#PBS -l mem=%dGB
#PBS -l walltime=%s
#PBS -l wd
#PBS -l storage=scratch/%s+gdata/%s

module load nci-parallel/1.0.0a
module load R/4.3.1

export ncores_per_task=1
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file cmds.txt --timeout 3600 --status status.txt --output-dir logs
', project, ncpus, mem_gb, walltime, project, project)

con2 <- file("job.sh", open = "wb")
writeLines(job_script, con2, sep = "\n")
close(con2)

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

cat(sprintf("%d tasks written to cmds.txt\n", length(cmds)))