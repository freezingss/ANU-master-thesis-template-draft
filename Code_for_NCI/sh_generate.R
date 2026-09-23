project <- "vk72"
ncpus <- 4
mem_gb <- 16
walltime <- "08:00:00"

JQ_grid <- list(c(50, 100), c(50, 150), c(50, 200), c(50, 250), c(50, 300),
                c(100, 100), c(100, 150), c(100, 200), c(100, 250), c(100, 300),
                c(150, 100), c(150, 150), c(150, 200), c(150, 250), c(150, 300),
                c(200, 100), c(200, 150), c(200, 200), c(200, 250), c(200, 300)
                )
seeds <- 1:10
methods <- c("wb")
K <- 2

cmds <- character(0)
for (jq in JQ_grid) {
  for (s in seeds) {
    for (m in methods) {
      cmds <- c(cmds, sprintf("Rscript run_one_cell.R %d %d %d %s %d", jq[1], jq[2], s, m, K))
    }
  }
}

cmds <- gsub("\r", "", cmds)

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
module load R/4.5.0

export ncores_per_task=1
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file cmds.txt --timeout 3600 --status status.txt --output-dir logs
', project, ncpus, mem_gb, walltime, project, project)

job_script <- gsub("\r", "", job_script)

con2 <- file("job.sh", open = "wb")
writeLines(job_script, con2, sep = "\n")
close(con2)

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

cat(sprintf("%d tasks written to cmds.txt\n", length(cmds)))