project <- "vk72"
ncpus <- 48
mem_gb <- 192
walltime <- "24:00:00"

JQ_grid <- expand.grid(
  J = c(100),
  Q = c(50, 100, 150, 200)
)

seeds <- 1:10
methods <- c("glmmTMB")
K <- 2

cmds <- character(0)
for (i in seq_len(nrow(JQ_grid))) {
  J_i <- JQ_grid$J[i]
  Q_i <- JQ_grid$Q[i]
  for (s in seeds) {
    for (m in methods) {
      cmds <- c(cmds, sprintf("Rscript run_one_cell_glmmTMB_test.R %d %d %d %s %d", J_i, Q_i, s, m, K))
    }
  }
}

cmds <- gsub("\r", "", cmds)

con <- file("cmds_glmmTMB_seed1_10.txt", open = "wb")
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
module load intel-compiler-llvm/2025.3.2
module load intel-mkl/2025.0.1
module load R/4.5.0

export ncores_per_task=4
export ncores_per_numanode=12

mpirun -np $((PBS_NCPUS/ncores_per_task)) --map-by ppr:$((ncores_per_numanode/ncores_per_task)):NUMA:PE=${ncores_per_task} nci-parallel --input-file cmds_glmmTMB_seed1_10.txt --timeout 14400 --status status_glmmTMB_seed1_10.txt --output-dir logs
', project, ncpus, mem_gb, walltime, project, project)

job_script <- gsub("\r", "", job_script)

con2 <- file("job_glmmTMB_seed1_10.sh", open = "wb")
writeLines(job_script, con2, sep = "\n")
close(con2)

dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

cat(sprintf("%d tasks written to cmds_glmmTMB_seed1_10.txt\n", length(cmds)))
