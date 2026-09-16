# Change the name of several function, make sure their effect is the same as before, and the cited name is correct  
# Notice: our own the version of R is R/4.5.1, does not match the required version of NCI (whose newest R version is R/4.5.0)

# # For VSC R language server
# install.packages("languageserver")
# install.packages("jsonlite")

# Used Packages
# install.packages("glmmTMB")
library(glmmTMB)
# install.packages("gllvm")
library(gllvm)
# install.packages("bench") # timing
library(bench)
# install.packages("nnet")
library(nnet)
# install.packages("Matrix") # Cholesky decomposition
library(Matrix)
# install.packages("RSpectra") # eigen decomposition
library(RSpectra)
# install.packages("ggplot2")
library(ggplot2)