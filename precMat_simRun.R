setwd("/Users/berichzinsou-daho/Documents/AMSE/PHD/Master thesis project/Financial Econometrics/Simulations")

####******** Estimation methods simulation********

source("sim_packages.R")
source("simulation_main.R")

stocks <- read.csv("stocks_prices.csv",sep = ",")
#stocks_return <- calculate_monthly_returns(stocks,600)

# Get factor from website and download it ans save as text file. then upload
file_lines <- read_lines("F-F_Research_Data_Factors.txt")
mkt_data <- read_table(file_lines)
mkt_fac <- data_frame(mkt_data)
mkt_fac$Date<- as.Date(paste0(mkt_fac$Date, "01"), format = "%Y%m%d")
mkt_fac$Date <- format(mkt_fac$Date, "%Y-%m")
mkt_fac <- mkt_fac %>% dplyr::select(c("Date","Mkt-RF","SMB","HML"))
#merge return and factor on date
#merged_data<- inner_join(stocks_return,mkt_fac, by = "Date")
#return_matrix <- merged_data %>% dplyr::select(-c("Date","Mkt-RF","SMB","HML","RF"))
#factor_matrix <- merged_data %>% dplyr::select(c("Mkt-RF","SMB","HML"))

## Simulate return. and factors
#sim_return <- simulateAssetReturns(return_matrix,factor_matrix,rho = 0.5, n = 400)
#save(sim_return,file ="sim_n400/simulatedReturn_n400_p600.RData")

##
config = read_yaml("sim_config.yaml")
N = config$number_of_sim # number of simulation
rho_1 = config$target_expected_return_MMV
#load("sim_n100/simulatedReturn_n100_p150.RData")

#y = sim_return$SimulatedReturns
#fac = sim_return$Factors
#covMat = sim_return$covMat
#TruePrecision = solve(covMat)

# Loop for 1000 replication
results_list <- list()
time_taken <- system.time({
  for (i in 1:N) {
    stocks_return <- calculate_monthly_returns(stocks,150)
    
    merged_data<- inner_join(stocks_return,mkt_fac, by = "Date")
    return_matrix <- merged_data %>% dplyr::select(-c("Date","Mkt-RF","SMB","HML"))
    factor_matrix <- merged_data %>% dplyr::select(c("Mkt-RF","SMB","HML"))
    
    
    sim_return <- simulateAssetReturns(return_matrix,factor_matrix,rho = 0.5, n = 100)
    
    y = sim_return$SimulatedReturns
    fac = sim_return$Factors
    covMat = sim_return$covMat
    TruePrecision = solve(covMat)
    
    results_list[[i]] <- AutoNodewisePrecMatrix (fac,y)
  }
})
# Save all matrices in a single .RData file
save(results_list, file = "BoopSim/sim_n100/sim_p150/AutoNodewise_prec_mat.RData")


load("BoopSim/sim_n100/sim_p150/AutoNodewise_prec_mat.RData")

SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- MaxOosSharpeRatio(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MosSR_AutoNodewise.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- SharpeRatioGMV(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "BoopSim/sim_n100/sim_p150/perfMeasure/GmvSR_AutoNodewise.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- SharpeRatioMMV(results_list[[i]],y,1/100)
  }
})
save(SharpeRatio_list, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MmvSR_AutoNodewise.RData")


SharpeRatio_list= list()
time_taken <- system.time({
  for (i in seq_len(N)) {
    SharpeRatio_list[[i]] <- MsrPortfolio(y,results_list[[i]])
  }
})
save(SharpeRatio_list, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MsrSR_AutoNodewise.RData")

# Evaluation by mean square error (mse)
#Evaluation for the performance measurement (sharpe ratio)
load("BoopSim/sim_n100/sim_p150/perfMeasure/MosSR_AutoNodewise.RData")
TrueSharpeRatio <- MaxOosSharpeRatio(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_MosSR_AutoNodewise.RData")

TrueSharpeRatio <- MaxOosSharpeRatio(y,TruePrecision)
mseSR <- numeric()
for (i in seq_len(N)) {
  mseSR <- {mean(apply(abs(TrueSharpeRatio- SharpeRatio_list[[i]]),1,sum))} 
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_abs_MosSR_AutoNodewise.RData")



load("BoopSim/sim_n100/sim_p150/perfMeasure/GmvSR_AutoNodewise.RData")
TrueSharpeRatio <- SharpeRatioGMV(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_GmvSR_AutoNodewise.RData")

TrueSharpeRatio <- SharpeRatioGMV(y,TruePrecision)
mseSR <- numeric()
for (i in seq_len(N)) {
  mseSR <- {mean(apply(abs(TrueSharpeRatio- SharpeRatio_list[[i]]),1,sum))} 
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_abs_GmvSR_AutoNodewise.RData")



load("BoopSim/sim_n100/sim_p150/perfMeasure/MmvSR_AutoNodewise.RData")
TrueSharpeRatio <- SharpeRatioMMV(TruePrecision,y,1/100)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_MmvSR_AutoNodewise.RData")

TrueSharpeRatio <- SharpeRatioMMV(TruePrecision,y,1/100)
mseSR <- numeric()
for (i in seq_len(N)) {
  mseSR <- {mean(apply(abs(TrueSharpeRatio- SharpeRatio_list[[i]]),1,sum))} 
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_abs_MmvSR_AutoNodewise.RData")



load("BoopSim/sim_n100/sim_p150/perfMeasure/MsrSR_AutoNodewise.RData")
TrueSharpeRatio <- MsrPortfolio(y,TruePrecision)
mseSR <- 0
for (i in seq_len(N)) {
  mseSR <- mseSR + (TrueSharpeRatio - SharpeRatio_list[[i]])^2
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_MsrSR_AutoNodewise.RData")

TrueSharpeRatio <- MsrPortfolio(y,TruePrecision)
mseSR <- numeric()
for (i in seq_len(N)) {
  mseSR <- {mean(apply(abs(TrueSharpeRatio- SharpeRatio_list[[i]]),1,sum))} 
}
mseSR <- mseSR / (N)
save(mseSR, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_abs_MsrSR_AutoNodewise.RData")


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
save(mse, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_AutoNodewise.RData")

mse_em <- 0
for (i in seq_len(N)) {
  mse_em <- mse_em + sum(TruePrecision - results_list[[i]])
}
mse_em <- mse_em / (n_true * m_true *N)
save(mse_em, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_EM_AutoNodewise.RData")

mse_ema <- 0
for (i in seq_len(N)) {
  mse_ema <- mse_ema + sum(abs(TruePrecision - results_list[[i]]))
}
mse_ema <- mse_ema / (n_true * m_true *N)
save(mse_ema, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_EMA_AutoNodewise.RData")

#norm<-function(A,B){max(abs(A-B))}
p = ncol(TruePrecision)
mse_norm <- numeric()
for (i in seq_len(N)) {
  mse_diff <- max(abs(TruePrecision-results_list[[i]])) 
  mse_norm[i] <- mse_diff / p^2
}
mse_norm_mean <- mean(mse_norm)
save(mse_norm, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_Norm_AutoNodewise.RData")
save(mse_norm_mean, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_Norm_mean_AutoNodewise.RData")

#norm_2<-function(A,B){mean(apply(abs(A-B),1,sum))}
p = ncol(TruePrecision)
mse_norm2 <- numeric()
for (i in seq_len(N)) {
  mse_diff2 <- {mean(apply(abs(TruePrecision-results_list[[i]]),1,sum))} 
  mse_norm2[i] <- mse_diff2 / p^2
}
mse_norm2_mean <- mean(mse_norm2)
save(mse_norm2, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_Norm2_AutoNodewise.RData")
save(mse_norm2_mean, file = "BoopSim/sim_n100/sim_p150/perfMeasure/MSE_Norm2_mean_AutoNodewise.RData")


