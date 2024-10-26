setwd("/Users/berichzinsou-daho/Documents/AMSE/PHD/Master thesis project/Financial Econometrics/Simulations")

####******** Estimation methods simulation********

source("sim_packages.R")
source("simulation_main.R")

##
config = read_yaml("sim_config.yaml")
N = config$number_of_sim # number of simulation
rho_1 = config$target_expected_return_MMV
load("sim_n100/simulatedReturn_n100_p50.RData")

y = sim_return$SimulatedReturns
fac = sim_return$Factors
covMat = sim_return$covMat
TruePrecision = solve(covMat)

# Loop for 1000 replication
results_list <- list()
time_taken <- system.time({
  for (i in 1:N) {
    
    results_list[[i]] <- ic.prec.matrix(fac,y,crit = "bic")
  }
})
# Save all matrices in a single .RData file
save(results_list, file = "sim_n100/sim_p150/NodewiseBic_prec_mat.RData")


load("sim_n100/sim_p150/NodewiseBic_prec_mat.RData")

SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- MaxOosSharpeRatio(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "sim_n100/sim_p150/perfMeasure/MosSR_NodewiseBic.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- SharpeRatioGMV(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "sim_n100/sim_p150/perfMeasure/GmvSR_NodewiseBic.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- SharpeRatioMMV(results_list[[i]],y,1/100)
  }
})
save(SharpeRatio_list, file = "sim_n100/sim_p150/perfMeasure/MmvSR_NodewiseBic.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- MsrPortfolio(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "sim_n100/sim_p150/perfMeasure/MsrSR_NodewiseBic.RData")

# Evaluation by mean square error (mse)
#Evaluation for the performance measurement ( sharpe ratio)
load("sim_n100/sim_p150/perfMeasure/MosSR_NodewiseBic.RData")
TrueSharpeRatio <- MaxOosSharpeRatio(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "sim_n100/sim_p150/perfMeasure/MSE_MosSR_NodewiseBic.RData")

load("sim_n100/sim_p150/perfMeasure/GmvSR_NodewiseBic.RData")
TrueSharpeRatio <- SharpeRatioGMV(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "sim_n100/sim_p150/perfMeasure/MSE_GmvSR_NodewiseBic.RData")

load("sim_n100/sim_p150/perfMeasure/MmvSR_NodewiseBic.RData")
TrueSharpeRatio <- SharpeRatioMMV(TruePrecision,y,1/100)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "sim_n100/sim_p150/perfMeasure/MSE_MmvSR_NodewiseBic.RData")

load("sim_n100/sim_p150/perfMeasure/MsrSR_NodewiseBic.RData")
TrueSharpeRatio <- MsrPortfolio(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "sim_n100/sim_p150/perfMeasure/MSE_MsrSR_NodewiseBic.RData")


# Evaluation for the precision matrix
# Dimensions of the matrices (assuming all matrices are the same size)
n_true <- nrow(TruePrecision)
m_true <- ncol(TruePrecision)

# Calculate MSE
mse <- 0
for (i in seq_len(N)) {
  mse <- mse + sum((TruePrecision - results_list[[i]])^2)
}
mse <- mse / (n_true * m_true *N)
save(mse, file = "sim_n100/sim_p150/perfMeasure/MSE_NodewiseBic.RData")
