setwd("/Users/berichzinsou-daho/Documents/AMSE/PHD/Master thesis project/Financial Econometrics/Simulations")
# Function to check and install packages
check_and_install_packages <- function(package_list) {
  for (package_name in package_list) {
    if (!require(package_name, character.only = TRUE)) {
      install.packages(package_name)
      library(package_name, character.only = TRUE)
    }
  }
}

packages_needed <- c("quantmod", "xts", "zoo", "rvest", "tidyverse", 
                     "tidyquant", "janitor", "purrr", "PerformanceAnalytics", 
                     "timetk", "MASS", "glmnet", "lars", "yaml")

# Define the function to get stock prices 
get_stocks_prices <- function(from_date, to_date, index_to_extract) {
                # index in quote
                # date in cote
  # Fetch index's stock symbols
  message("Fetching stocks symbols...")  # message to indicate that stocks symbols are fetching
  stocks <- tq_index(index_to_extract)
  stock_sym <- as.character(stocks$symbol)
  
  # Use map_df from purrr to iterate over symbols and fetch prices
  # This approach is more efficient and concise than a loop
  message("Fetching stock prices for each symbol...")     #message to indicate that stocks prices are fetching
  stocks_price <- purrr::map_df(stock_sym, function(symbol) {
    tryCatch({
      tq_get(symbol,
             get = "stock.prices",
             from = from_date,
             to = to_date) %>%
        dplyr::mutate(Symbol = symbol) # Add a symbol column for identification
    }, error = function(e) {
      NULL # In case of error, return NULL 
    })
  }, .id = "Symbol")
  
  return(stocks_price)
}


# Get factor from website and download it ans save as text file. then upload
file_lines <- read_lines("F-F_Research_Data_Factors.txt")
mkt_data <- read_table(file_lines)
mkt_fac <- data_frame(mkt_data)
mkt_fac$Date<- as.Date(paste0(mkt_fac$Date, "01"), format = "%Y%m%d")
mkt_fac$Date <- format(mkt_fac$Date, "%Y-%m")


# Parameter to be set 
config = read_yaml("sim_config.yaml")
n = config$ number_of_obs  # number of observation
p = config$number_of_asset   # number of asset

# we select first p asset from those we've get

# Define the function
calculate_monthly_returns <- function(x,p) {
                            # x : data on which we would like to work
                            #p = is the number of asset "p " we would like to pick for the simulation
  
    # Filter stock prices data for selected symbols
     
     adjusted <- x[,c("date","symbol","adjusted")]
     adjusted["date"] = as.Date(adjusted[["date"]])
     stcks <- adjusted %>% 
          pivot_wider(
          names_from = symbol, 
          values_from = adjusted
        )
     new_df <- stcks %>% select_if(~ !any(is.na(.)))
     
     # Randomly select p symbols
     selected_symbols <- sample(colnames(new_df[,-1]), size = p, replace = TRUE)
     columns_to_select <- c("date", selected_symbols)
     #a <- new_df %>% dplyr::select("date",all_of(selected_symbols))
     #a <- cbind(new_df["date"],new_df[,selected_symbols]) 
     a <- new_df[, columns_to_select]
     #a <- as.xts(a)
     prices_monthly <- to.monthly(a, indexAt = "last", OHLC = FALSE)

     prices_monthly <- prices_monthly %>% 
       mutate(across(everything(), ~as.numeric(trimws(.))))
  
     asset_returns_xts <- prices_monthly %>% Return.calculate( method = "log") %>% na.omit()
     
     #asset_returns_xts <- na.omit(Return.calculate(prices_monthly, method = "log"))
     returns <- rownames_to_column(asset_returns_xts)
     returns <- rename(returns, "Date" = "rowname")
     returns$Date <- as.Date(returns$Date)
     ret <- returns
     ret$Date <- format(ret$Date,"%Y-%m")
  return(ret)
}


simulateAssetReturns <- function(assets_returns, ret_factor, rho, n) {
            #assets_returns : it's about the matrix of assets returns only, without date columns
            #ret_factor : here need matrix or vector of factor
            #rho # the toeplitz coefficient that we need to compute th matrx for the covariance structure
            #n : number of observation we would like to generate
  
  # Assuming sim_data[sel_col] and FF_factor are defined outside and passed as parameters
  X <- as.matrix(ret_factor) # Fama & French factor as passed parameter
  y <- as.matrix(assets_returns)
  
  # Run linear model
  dgp_1 <- lm(y ~ X)
  betas_j <- dgp_1$coefficients[-1,]
  
  # Generate idiosyncratic error
  e <- dgp_1$residuals
  sigma_e <- cov(e)
  tptz_mtx <- toeplitz((rho)^(0:(ncol(e)-1))) # generate Toeplitz matrix
  e_hat <- sigma_e * tptz_mtx # covariance matrix of error
  mean_vector <- rep(0, ncol(e))
  err <- mvrnorm(n, mu = mean_vector, Sigma = e_hat) # Generate errors
  
  # Generate factors
  f_mean_vector <- colMeans(ret_factor) # factor mean vector
  f_1 <- mvrnorm(n, mu = f_mean_vector, Sigma = cov(ret_factor))
  
  # Generate simulated data
  X <- as.matrix(f_1)  # Matrix of factors
  Y <- as.matrix(X) %*% betas_j + err
  #Y <- as.matrix(cbind(1, X)) %*% betas_j + err
  return(list(Factors = X, SimulatedReturns = Y,covMat = e_hat))
}
##
#
###
ic.prec.matrix = function (x, y, crit=c("bic","aic","aicc","hqc","gic"),alpha = 1)
{
  # y : is the matrix of return
  # x : is the matrix of factor
  # if there is no factor , then create a matrix of zero in order to run the function
  #if (!is.matrix(x)) x <- as.matrix(x)
  #if (!is.matrix(y)) stop("y needs to be a matrix.")
  x = as.matrix(x)
  y = as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  omega_hat_gic <- matrix(nrow = p, ncol = p)
  
  for (col_index in 1 : ncol(y)) 
    
  {
  
  # Select a column from the matrix
  y_j <- as.matrix(y[,col_index])  # return of the jth asset
  Y_minusj  = as.matrix(y[,-col_index])  ## return of all asset without the jth
  
  ### start nodewise and precision matrix
  step_1 <- lm(y_j ~ x) # regress the jth asset on the generated factor
  û_j <- as.matrix(step_1$residuals)
  
  step_2 <- lm(Y_minusj ~ x)
  U_minusj <- as.matrix(step_2$residuals)
  
  choosen_crit=match.arg(crit)
  n=nrow(y)
  model = glmnet(U_minusj , û_j ,alpha = alpha)
  coef = coef(model)
  lambda = model$lambda
  df = model$df
  
  if(alpha==0){
    xs = scale(x)
    I = diag(ncol(x))
    xx = t(xs)%*%xs
    for(i in 1:length(lambda)){
      aux = solve(xx + I*lambda[i])
      df[i] = sum(diag(xs%*%aux%*%t(xs)))
    }
    
  }
  
  û_jhat=cbind(1,U_minusj)%*%coef
  residuals = (û_j - û_jhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)
  
  nvar = df + 1
  bic = n*log(mse)+nvar*log(n)
  aic = n*log(mse)+2*nvar
  aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
  hqc = n*log(mse)+2*nvar*log(log(n))
  gic = (sse/n) + nvar*(log(p-1)*(log(log(n))/n))
  
  crit_use=switch(choosen_crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc,gic=gic)
  
  selected_mod=best.model = which(crit_use == min(crit_use))
  
  ic=c(bic=bic[selected_mod],aic=aic[selected_mod],aicc=aicc[selected_mod],hqc=hqc[selected_mod],gic=gic[selected_mod])
  
  #precision matrix with information criterium lambda
  gamma_j_gic <- as.matrix(coef[,selected_mod][-1]) 
  tho_2j_gic <- t(û_j)%*%(û_j - U_minusj%*%gamma_j_gic)/n
  
  c_j_gic <- matrix(0, nrow = 1, ncol = ncol(y))
  c_j_gic[,col_index] <- 1
  c_j_gic[c_j_gic == 0] <- t(-gamma_j_gic)
  omega_hat_j_gic <- 1/as.numeric(tho_2j_gic) * c_j_gic
  omega_hat_gic[col_index,] <-omega_hat_j_gic
  
  
  }
  omega_hat_gic <- (omega_hat_gic+t(omega_hat_gic))/2
  return(omega_hat_gic)
}
###############

cv.prec.matrix = function (x, y,lambda_crit=c("lambda.min","lambda.1se"))
{
  # y : is the matrix of return
  # x : is the matrix of factor
  # if there is no factor , then create a matrix of zero in order to run the function
  #if (!is.matrix(x)) x <- as.matrix(x)
  #if (!is.matrix(y)) stop("y needs to be a matrix.")
  x = as.matrix(x)
  y = as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  omega_hat_cv <- matrix(nrow = p, ncol = p)
  
  for (col_index in 1 : ncol(y)) 
    
  {
    
    # Select a column from the matrix
    y_j <- as.matrix(y[,col_index])  # return of the jth asset
    Y_minusj  = as.matrix(y[,-col_index])  ## return of all asset without the jth
    
    ### start nodewise and precision matrix
    step_1 <- lm(y_j ~ x) # regress the jth asset on the generated factor
    û_j <- as.matrix(step_1$residuals)
    
    step_2 <- lm(Y_minusj ~ x)
    U_minusj <- as.matrix(step_2$residuals)
    
    
    chosen_crit <- match.arg(lambda_crit)
    #
    #perform k-fold cross-validation to find optimal lambda value
    cv_4_lamda <- cv.glmnet(U_minusj,û_j, alpha = 1) # cross validation with K-fole by default
    lambda_used <- if (chosen_crit == "lambda.min") cv_4_lamda$lambda.min else cv_4_lamda$lambda.1se
    
    
    #lambda_j <- cv_4_lamda$lambda.min # get the lambda that minimise MSE
    ## nodewise and precision matrix with cv lambda
    nd_wise <- glmnet(U_minusj,û_j,alpha = 1, lambda = lambda_used) #nodewise estimator for asset j
    gamma_j <- as.matrix(nd_wise$beta) # estimate coef from the nodewise
    tho_2j <- t(û_j)%*%(û_j - U_minusj%*%gamma_j)/n
    
    c_j <- matrix(0, nrow = 1, ncol = p)
    c_j[,col_index] <- 1
    c_j[c_j == 0] <- t(-gamma_j)
    omega_hat_j <- 1/as.numeric(tho_2j) * c_j
    omega_hat_cv[col_index,] <-omega_hat_j
   
  }
  omega_hat_cv <- (omega_hat_cv+t(omega_hat_cv))/2
  return(omega_hat_cv)
}

#
calculate_gamma_hat <- function(omega_hat, x, y) {
  # gamma hat here stand for the precision matrix of all asset in the portfolio
   # y : matrix of return
   # x : matrix of factor
   # omega_hat : precision matrix of error
  
  omega_hat_sym <- (omega_hat+t(omega_hat))/2
  B_hat_reg <- lm(y ~ x) # estimation for all asset return
  B_hat <- B_hat_reg$coefficients # coeficients foe each association return-factor
  
  # Ensure inputs are matrices
  omega_hat <- as.matrix(omega_hat)
  B_hat <- as.matrix(B_hat)
  cov_ft <- cov(x)
  omega_hat_sym <- as.matrix(omega_hat_sym)
  
  # Part 1
  part1 <- t(B_hat[-1, ]) %*% omega_hat
  
  # Part 2
  part2 <- solve(cov_ft) + B_hat[-1, ] %*% omega_hat_sym %*% t(B_hat[-1, ])
  
  # Combining and completing the operation
  gamma_hat <- omega_hat - (t(part1) %*% (solve(part2) %*% t(omega_hat %*% t(B_hat[-1, ]))))
  
  return(gamma_hat)
}

#Right formula to check results with, if necessary
#gamma_hat <- omega_hat - (t(B_hat[-1,]%*%omega_hat) %*% (solve(solve(cov_ft) +
                                                                # B_hat[-1,]%*%omega_hat_sym%*%t(B_hat[-1,]))) %*% t(omega_hat%*%t(B_hat[-1,])))


##******Precision matrix with POET*******

calculatePoetGammaHat <- function(y) {
  # Calculate the optimal number of factors
  opt_k <- POETKhat(y)
  
  # Calculate the optimal threshold
  opt_thr <- POETCmin(y, opt_k, method = 'soft', estimator = 'vad')
  
  # Calculate the covariance matrix using POET
  poet_cov <- POET(t(y), opt_k, opt_thr, method = 'soft', estimator = 'vad')$SigmaY
  
  # Calculate the precision matrix (inverse of the covariance matrix)
  poet_gamma_hat <- solve(poet_cov)
  
  # Return the rounded precision matrix
  return(round(poet_gamma_hat, 5))
}

##******Precision matrix with non linear shrinkage####

##
#function to check against the formula at the first use
maxser_weigth <- function(y) {
  # Ensure Y is a matrix, mu_hat is a matrix or vector as appropriate
  y <- as.matrix(y)
  n = nrow(y)
  p = ncol(y)
  mu_hat =  matrix(nrow = 1 , ncol = ncol(y))
  for (i in 1:ncol(y)) {
    mu_hat[i] <- mean(y[,i ])
  }
  one_vec_p <- matrix(1, nrow = ncol(y), ncol = 1)
  mu_hat <- as.matrix(mu_hat)
  
  # Calculate sigma_sample as the inverse of the covariance matrix of Y
  sigma_sample <- solve(cov(y))
  
  # Compute theta_sample using mu_hat, sigma_sample
  theta_sample <- mu_hat %*% sigma_sample %*% t(mu_hat)
  
  # Calculate theta_hat_mser based on the formula provided
  theta_hat_mser <- (((n - p - 2) * theta_sample) - p) / n
  
  r_c <- (1 + theta_hat_mser) / sqrt(as.numeric(theta_hat_mser))#return maxser method based for unsconstrained optimization
  
  # Create r_c_vec as a matrix
  r_c_vec <- matrix(r_c, ncol = 1, nrow = n)
  
  # Perform penalized regression (use lars to perform the lasso) to get the weights
  lars_model <- lars(y, r_c_vec, type = 'lasso')
  omega_mser_hat <- lars_model$beta # Extract the estimated weights
  
  return(omega_mser_hat)
}

###
SharpeRatioGMV <- function(y, gamma_hat) {
  # gamma hat : precision matrix of error or for asset depending on whether we use factor or not
  # Ensure Y and gamma_hat are matrices
  y <- as.matrix(y)
  gamma_hat <- as.matrix(gamma_hat)
  
  p <- ncol(y)  # Number of assets
  
  # Calculate mean returns for each asset
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:p) {
    mu_hat[i] <- mean(y[, i])
  }
  
  # Create a vector of ones with length equal to the number of assets
  one_vec_p <- matrix(1, nrow = ncol(y), ncol = 1)
  
  # Calculate the Sharpe ratio under GMV
  gmv_sharpeRatio <- sqrt(p) * (t(one_vec_p) %*% gamma_hat %*% t(mu_hat) / p) * 
    as.numeric(((t(one_vec_p) %*% gamma_hat %*% one_vec_p) / p)^(-1/2))
  
  return(gmv_sharpeRatio)
}

GMVweights <- function(y, gamma_hat) {
  # gamma hat : precision matrix of error or for asset depending on whether we use factor or not
  # Ensure Y and gamma_hat are matrices
  y <- as.matrix(y)
  gamma_hat <- as.matrix(gamma_hat)
  
  p <- ncol(y)  # Number of assets
  
  # Calculate mean returns for each asset
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:p) {
    mu_hat[i] <- mean(y[, i])
  }
  
  # Create a vector of ones with length equal to the number of assets
  one_vec_p <- matrix(1, nrow = ncol(gamma_hat), ncol = 1)
  
  # Calculate the Sharpe ratio under GMV
  gmv_weigths <-  (gamma_hat %*% one_vec_p) / as.numeric(t(one_vec_p) %*% gamma_hat %*% one_vec_p)
  
  return(gmv_weigths)
}
##

SharpeRatioMMV <- function(gamma_hat,y,rho) {
  #rho : target expected return
  p = ncol(y)
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:ncol(y)) {
    mu_hat[i] <- mean(y[, i])
  }
  # Ensure gamma_hat and mu_hat are matrices/vectors
  gamma_hat <- as.matrix(gamma_hat)
  mu_hat <- as.matrix(mu_hat)
  
  # Ensure the dimensions and orientations are correct
  #if (ncol(mu_hat) != 1) {
   # mu_hat <- t(mu_hat)  # Transpose mu_hat if it's not a column vector
  #}
  
  # Create a vector of ones (p x 1)
  one_vec_p <- matrix(1, nrow = p, ncol = 1)
  
  # Calculate A_hat, F_hat, and D_hat
  A_hat <- as.numeric((t(one_vec_p) %*% gamma_hat %*% one_vec_p) / p)
  F_hat <- (t(one_vec_p) %*% gamma_hat %*% t(mu_hat)) / p
  D_hat <- (mu_hat %*% gamma_hat %*% t(mu_hat)) / p
  
  # Calculate the Sharpe ratio under MMV
  MMV_sharpeRatio <- rho * (sqrt(p * (A_hat * D_hat - F_hat^2) / (A_hat * rho^2 - 2 * F_hat * rho + D_hat)))
  
  return(MMV_sharpeRatio)
}
#
##MMV weights
MMVweight <- function(gamma_hat,y,rho) {
  #rho : target expected return
  p = ncol(y)
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:ncol(y)) {
    mu_hat[i] <- mean(y[, i])
  }
  # Ensure gamma_hat and mu_hat are matrices/vectors
  gamma_hat <- as.matrix(gamma_hat)
  mu_hat <- as.matrix(mu_hat)
  # Create a vector of ones (p x 1)
  one_vec_p <- matrix(1, nrow = p, ncol = 1)
  
  # Calculate A_hat, F_hat, and D_hat
  A_hat <- as.numeric((t(one_vec_p) %*% gamma_hat %*% one_vec_p) / p)
  F_hat <- (t(one_vec_p) %*% gamma_hat %*% t(mu_hat)) / p
  D_hat <- (mu_hat %*% gamma_hat %*% t(mu_hat)) / p
  
  # Calculate the Sharpe ratio under MMV
  MMV_weight <- as.numeric(((D_hat-rho*F_hat)/(A_hat*D_hat-F_hat^2))) * (t(one_vec_p) %*% gamma_hat/p) +
   as.numeric (((rho*A_hat - F_hat)/(A_hat*D_hat - F_hat^2))) * (mu_hat %*% gamma_hat /p)
 
  return(MMV_weight)
}

###

# function to check against the formula at the first use
MaxOosSharpeRatio <- function(y, gamma_hat) {
  # gamma hat : precise matrix of assets or for error if there is not factor
  # y : return
  p = ncol(y)
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:p) {
    mu_hat[i] <- mean(y[, i])
  }
  # Ensure inputs are in the correct matrix/vector format
  y <- as.matrix(y)
  mu_hat <- as.matrix(mu_hat)
  gamma_hat <- as.matrix(gamma_hat)
  
  # Calculate the numerator (portfolio's expected return squared)
  expected_return_squared <- mu_hat %*% gamma_hat %*% t(mu_hat)
  
  # Calculate the denominator (standard deviation of portfolio's return)
  # Here we use the covariance matrix of Y, adjusted by gamma_hat
  variance <- mu_hat %*% t(gamma_hat) %*% cov(y) %*% gamma_hat %*% t(mu_hat)
  standard_deviation <- sqrt(variance)
  
  # Calculate the maximum out-of-sample Sharpe ratio
  MOS_sharpeRatio <- expected_return_squared / standard_deviation
  
  return(MOS_sharpeRatio)
}

#
# function to check against the formula at the first use
MsrPortfolio <- function(y, gamma_hat) {
  #function for maximum sharpe ratio portfolio
  # gamma hat : precise matrix of assets or for error if there is not factor
  # y : return
  p = ncol(y)
  mu_hat <- matrix(nrow = 1, ncol = p)
  for (i in 1:p) {
    mu_hat[i] <- mean(y[, i])
  }
  
  one_vec_p <- matrix(1, nrow = ncol(y), ncol = 1)
  # Ensure inputs are correctly formatted as matrices
  mu_hat <- as.matrix(mu_hat)
  gamma_hat <- as.matrix(gamma_hat)
  one_vec_p <- as.matrix(one_vec_p)
  
  # Calculate condition
  cdtion <- t(one_vec_p) %*% gamma_hat %*% t(mu_hat)
  
  # Calculate unconstrained Maximum Sharpe Ratio (MSR)
  Msr_nc <- mu_hat %*% gamma_hat %*% t(mu_hat)
  
  # Calculate constrained Maximum Sharpe Ratio (MSR)
  Msr_c <- mu_hat %*% gamma_hat %*% t(mu_hat) - 
    ((t(one_vec_p) %*% gamma_hat %*% t(mu_hat))^2) / (t(one_vec_p) %*% gamma_hat %*% one_vec_p)
  
  # Determine the optimal Maximum Sharpe Ratio based on condition
  MaxSharperatio_p <- ifelse(cdtion > 0, Msr_nc, Msr_c)
  
  return(MaxSharperatio_p)
}
#
##
CholeskyPrecisionMatrix <- function(y,lambda_crit=c("lambda.min","lambda.1se")) {
  # y : matrix of return
  n <- nrow(y)
  p <- ncol(y)
  
  # Initialize B and G
  b <- matrix(0, nrow = p, ncol = p)  # Coefficients matrix
  G <- diag(nrow = p)  # Variance of residuals, diagonal matrix for simplicity
  
  
  for (i in 1:p) {
    if (i == 1) {
      residuals <- y[, i] - mean(y[, i])
    } else if(i==2) {predictors <- y[, 1]
    response <- y[, i]
    
    
    model <- lm(response ~ predictors - 1) # '-1' to exclude the intercept
    
    # Store coefficients in b
    b[i, 1:(i-1)] <- coef(model)  # Coefficients are stored directly
    
    residuals <- model$residuals # compute residuals
   
    chosen_crit <- match.arg(lambda_crit)
    } else {
      predictors <- y[, 1:(i-1)]
      response <- y[, i]
      
      
      cvLambda <- cv.glmnet(as.matrix(predictors), response, alpha = 1) 
      
      # Use chosen criterion for lambda
      lambda_used <- if (chosen_crit == "lambda.min") cvLambda$lambda.min else cvLambda$lambda.1se
      
      lassoModel <- glmnet(as.matrix(predictors), response, alpha = 1, lambda = lambda_used)
      
      coef_matrix <- coef(lassoModel, s = lambda_used)
      b[i, 1:(i-1)] <- coef_matrix[-1, 1]  # Exclude the intercept
      
      predictions <- predict(lassoModel, newx = as.matrix(predictors), s = lambda_used)
      residuals <- response - predictions
    }
    
    G[i, i] <- var(residuals)
  }
  
  I <- diag(p)
  B <- I - b
  G = 1 / diag(G)
  chlesky_precMatrix <- t(B) %*% diag(G) %*% B
  
  #return(list(b=b,B = B, G = G, chlesky_precMatrix = chlesky_precMatrix))
  return(chlesky_precMatrix = chlesky_precMatrix)
}
### cholesky precision matrix with lambda selection:
##
ic.CholeskyPrecisionMatrix <- function(y, crit=c("bic","aic","aicc","hqc","gic")) {
  # y : matrix of return
  n <- nrow(y)
  p <- ncol(y)
  
  # Initialize B and G
  b <- matrix(0, nrow = p, ncol = p)  # Coefficients matrix
  G <- diag(nrow = p)  # Variance of residuals, diagonal matrix for simplicity
  
  
  for (i in 1:p) {
    if (i == 1) {
      residuals <- y[, i] - mean(y[, i])
    } else if(i==2) {predictors <- y[, 1]
    response <- y[, i]
    
    
    model <- lm(response ~ predictors - 1) # '-1' to exclude the intercept
    
    # Store coefficients in b
    b[i, 1:(i-1)] <- coef(model)  # Coefficients are stored directly
    
    residuals <- model$residuals # compute residuals
    
    } else {
      
      choosen_crit=match.arg(crit)
      predictors <- y[, 1:(i-1)]
      response <- y[, i]
      
      lassoModel <- glmnet(as.matrix(predictors), response, alpha = 1)
      coef = coef(lassoModel)
      df = lassoModel$df
      #predictions=cbind(1,predictors)%*%coef
      predictions <- predict(lassoModel, newx = as.matrix(predictors))
      residual <- response - predictions
      
      mse = colMeans(residual^2)
      sse = colSums(residual^2)
      
      nvar = df + 1
      bic = n*log(mse)+nvar*log(n)
      aic = n*log(mse)+2*nvar
      aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
      hqc = n*log(mse)+2*nvar*log(log(n))
      gic = (sse/n) + nvar*(log(p-1)*(log(log(n))/n))
      
      crit_use=switch(choosen_crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc,gic=gic)
      
      selected_mod=best.lassoModel = which(crit_use == min(crit_use))
      
      ic=c(bic=bic[selected_mod],aic=aic[selected_mod],aicc=aicc[selected_mod],hqc=hqc[selected_mod],gic=gic[selected_mod])
      
      coef_matrix <- as.matrix(coef[,selected_mod][-1])
      b[i, 1:(i-1)] <- coef_matrix  # Exclude the intercept
      
      coef1 <-  as.matrix(coef[,selected_mod]) 
      prediction=cbind(1,predictors)%*%coef1
      residuals <- response - prediction
    }
    G[i, i] <- var(residuals)
  }
  
  I <- diag(p)
  B <- I - b
  G = 1 / diag(G)
  chlesky_precMatrix <- t(B) %*% diag(G) %*% B
  
  #return(list(b=b,B = B, G = G, chlesky_precMatrix = chlesky_precMatrix))
  return(chlesky_precMatrix = chlesky_precMatrix)
}

nlshrinkPrecMatrix <- function(x) {
  
  # Perform the non-linear shrinkage estimation
  nl_LW <- nlshrink_cov(x, k = 0, method = "nlminb")
  
  # Solve the matrix (invert it)
  nl_LW_gamma_hat <- solve(nl_LW)
  nlshrink_precMatrix = nl_LW_gamma_hat
  
  return(nlshrink_precMatrix)
}


lshrinkPrecMatrix <- function(x) {
  
  # Perform the non-linear shrinkage estimation
  l_LW <- linshrink_cov(x, k = 0)
  
  #get the inverse
  l_LW_gamma_hat <- solve(l_LW)
  lshrink_precMatrix =  l_LW_gamma_hat

  return(lshrink_precMatrix)
}


source("Autometrics_Function.R")
Auto_CholeskyPrecisionMatrix <- function(y, alpha = 0.05) {
  # y : matrix of return
  n <- nrow(y)
  p <- ncol(y)
  
  # Initialize B and G
  b <- matrix(0, nrow = p, ncol = p)  # Coefficients matrix
  G <- diag(nrow = p)  # Variance of residuals, diagonal matrix for simplicity
  
  for (i in 1:p) {
    if (i == 1) {
      residuals <- y[, i] - mean(y[, i])
    } else {
      predictors <- y[, 1:(i-1)]
      response <- y[, i]
      
      if (i == 2) {
        model <- lm(response ~ predictors - 1) # '-1' excludes the intercept
        b[i, 1:(i-1)] <- coef(model)  # Store coefficients
        residuals <- model$residuals
      } else {
        # Autometrics
        Auto <- run_Auto_Ox(response, as.matrix(predictors), alpha = alpha)
        if (is.null(Auto$coefs)) {
          coef_matrix <- matrix(0, nrow = 1, ncol = i-1)
        } else {
          coef_matrix <- as.matrix(Auto$coefs)
        }
        
        if (!is.null(Auto$idx_var_sel) && length(Auto$idx_var_sel) > 0) {
          b[i, Auto$idx_var_sel] <- coef_matrix
          model <- lm(response ~ predictors[, Auto$idx_var_sel] - 1)
          residuals <- model$residuals
        }
      }
    }
    G[i, i] <- var(residuals)  # Update the variance of residuals
  }
  
  I <- diag(p)
  B <- I - b
  G <- diag(1 / diag(G))
  chlesky_precMatrix <- t(B) %*% G %*% B
  
  return(chlesky_precMatrix = chlesky_precMatrix)
}


AutoNodewisePrecMatrix = function (x,y,alpha = 0.05)
{
  # y : is the matrix of return
  # x : is the matrix of factor
  # if there is no factor , then create a matrix of zero in order to run the function
  #if (!is.matrix(x)) x <- as.matrix(x)
  #if (!is.matrix(y)) stop("y needs to be a matrix.")
  x = as.matrix(x)
  y = as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  omega_hat_cv <- matrix(nrow = p, ncol = p)
  for (col_index in 1 : ncol(y)) 
    
  {
    
    # Select a column from the matrix
    y_j <- as.matrix(y[,col_index])  # return of the jth asset
    Y_minusj  = as.matrix(y[,-col_index])  ## return of all asset without the jth
    
    ### start nodewise and precision matrix
    step_1 <- lm(y_j ~ x) # regress the jth asset on the generated factor
    û_j <- as.matrix(step_1$residuals)
    
    step_2 <- lm(Y_minusj ~ x)
    U_minusj <- as.matrix(step_2$residuals)
    
    # Autometrics
    Auto <- run_Auto_Ox(û_j, as.matrix(U_minusj), alpha = alpha)
    if (is.null(Auto$coefs)) {
      coef_matrix <- matrix(0, nrow = 1, ncol = ncol(U_minusj))  # Assuming `U_minusj` columns match expected
    } else {
      coef_matrix <- as.matrix(Auto$coefs)
    }
    gamma_j<- matrix(0, nrow = 1, ncol = p-1) 
    if (!is.null(Auto$idx_var_sel) && length(Auto$idx_var_sel) > 0) {
      gamma_j[, c(Auto$idx_var_sel)] <- coef_matrix
    }
    #gamma_j[, c(Auto$idx_var_sel)] <- coef_matrix
    tho_2j <- t(û_j)%*%(û_j - U_minusj%*%t(gamma_j))/n
    
    c_j <- matrix(0, nrow = 1, ncol = p)
    c_j[,col_index] <- 1
    c_j[c_j == 0] <- t(-gamma_j)
    omega_hat_j <- 1/as.numeric(tho_2j) * c_j
    omega_hat_cv[col_index,] <-omega_hat_j
    
  }
  return(omega_hat_cv)
}


