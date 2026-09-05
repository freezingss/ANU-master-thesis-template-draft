# Results Storage

# RESULTS_FILE <- "test_v2_freezed+eigen_results_v1.rds" from test_v2_freezed+eigen.R
# compare three models: base/lambda_correction/combo (lambda_correction + freezed sigma2) and give the output
results1 <- readRDS("test_v2_freezed+eigen_results_v1.rds")
summary(results1)