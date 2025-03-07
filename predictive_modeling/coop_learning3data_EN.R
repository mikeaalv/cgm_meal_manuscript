# model for SSPG
# Load necessary libraries
library(caret)
library(glmnet)

# Source any required custom functions
#source("./cooperative_regression_function.R")
source("cooperative_learning_multi_view.R")


# Load your new dataset
new_data <- read.csv('./data/predictive_data_ver4.csv')
required_columns <- c("sspg_avg_all", "insulin_fasting_avg_all", "bmi_avg_all","a1c_avg_all","fbg_avg_all","age_today_avg_all","sex",grep("^meta_", colnames(new_data), value = TRUE), grep("^lip_", colnames(new_data), value = TRUE))

new_data <- new_data[complete.cases(new_data[, required_columns]), ]


# Define the MSE calculation function
calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

# Assign your new data variables
y <- new_data$sspg_avg_all  # Outcome variable
x1 <- as.matrix(new_data[, c("insulin_fasting_avg_all", "bmi_avg_all","a1c_avg_all","fbg_avg_all","age_today_avg_all","sex")])
x2 <- as.matrix(new_data[, grep("^meta_", colnames(new_data))])
x3 <- as.matrix(new_data[, grep("^lip_", colnames(new_data))])


# Define x_fil and z_fil for filtered feature matrices
x1_fil <- x1
x2_fil <- x2
x3_fil <- x3

# Set model parameters
fit_mode = "min"
simN = 100
nfolds = 5
train_frac = 0.95
sim_seed = 123
val_frac = 0.4
set.seed(sim_seed)
alpha_value <-.9

alphalist = c(0.4,0.5,0.6,0.7,0.8,0.9,1)



# Initialize matrices to store errors and supports
err_train_coop = err_test_coop = support_coop = 
  diff_err = matrix(NA, simN, length(alphalist))
new_train_coop = new_cv_coop = new_test_coop = new_support_coop = matrix(NA, simN, length(alphalist))
new_train_coop_no_pf = new_support_coop_no_pf = new_test_coop_no_pf = new_cv_coop_no_pf = matrix(NA, simN, length(alphalist))

new_err_test_coop_validation = new_coop_support_validation = new_alpha_select = rep(NA, simN)

err_train_lasso = err_test_lasso = lasso_support = 
  X1_err_train_lasso = X1_err_test_lasso = X1_lasso_support = 
  X2_err_train_lasso = X2_err_test_lasso = X2_lasso_support = 
  X3_err_train_lasso = X3_err_test_lasso = X3_lasso_support = rep(NA, simN)


coop_selected_by_cv = support_by_cv = alpha_by_cv = rep(NA, simN)
coop_selected_by_cv_no_pf = support_by_cv_no_pf = alpha_by_cv_no_pf = rep(NA, simN)
err_null = rep(NA, simN)
err_fuse = rep(NA, simN)
support_fuse_late = rep(NA, simN)
coop_coeffs = rep(NA, simN)
coop_coeffs_no_pf = rep(NA, simN)

nonzero_coeffs <- data.frame(iteration=integer(),
                             model=character(),
                             alpha=numeric(),
                             coefficient_value=numeric(),
                             coefficient_name=character(),
                             stringsAsFactors=FALSE)
X1_test_pearson = X2_test_pearson = X3_test_pearson = 
  fuse_test_pearson = early_fusion_pearson = rep(NA, simN)

coop_no_pf_pearson_best = coop_pearson_best = rep(NA, simN)


# Simulation loop
for (ii in 1:simN){
  cat(ii, "\n")
  
  smp_size_train = floor(train_frac * nrow(x1_fil)) 
  train_ind = sort(sample(seq_len(nrow(x1_fil)), size = smp_size_train))
  test_ind = setdiff(seq_len(nrow(x1_fil)), train_ind)
  
  train_X1_raw <- x1_fil[train_ind, ]
  test_X1_raw <- x1_fil[test_ind, ]
  
  train_X2_raw <- x2_fil[train_ind, ]
  test_X2_raw <- x2_fil[test_ind, ]
  
  train_X3_raw <- x3_fil[train_ind, ]
  test_X3_raw <- x3_fil[test_ind, ]
  
  # Print dimensions for debugging
  cat("Initial dimensions of train_X1_raw:", dim(train_X1_raw), "\n")
  cat("Initial dimensions of train_X2_raw:", dim(train_X2_raw), "\n")
  cat("Initial dimensions of train_X3_raw:", dim(train_X3_raw), "\n")
  
  
  # Preprocess train_X1_raw
  preprocess_values_train_X1 = preProcess(train_X1_raw, method = c("center", "scale"))
  train_X1 = predict(preprocess_values_train_X1, train_X1_raw)
  test_X1 = predict(preprocess_values_train_X1, test_X1_raw)
  
  # Preprocess train_X2_raw
  preprocess_values_train_X2 = preProcess(train_X2_raw, method = c("center", "scale"))
  train_X2 = predict(preprocess_values_train_X2, train_X2_raw)
  test_X2 = predict(preprocess_values_train_X2, test_X2_raw)
  
  # Preprocess train_X3_raw
  preprocess_values_train_X3 = preProcess(train_X3_raw, method = c("center", "scale"))
  train_X3 = predict(preprocess_values_train_X3, train_X3_raw)
  test_X3 = predict(preprocess_values_train_X3, test_X3_raw)
  
  
  train_y <- y[train_ind]
  test_y <- y[test_ind]
  
  foldid = sample(rep_len(1:nfolds, dim(train_X1)[1]))
  
  # Null model
  print("Null model")
  err_null[ii] <- calc_mse(mean(train_y), test_y)
  print(err_null[ii])
  
  # Separate models
  # Separate models for X1
  print("Only X1")
  X1_lasso_fit = cv.glmnet(train_X1, train_y, alpha = alpha_value, standardize = F, foldid=foldid)
  
  if (fit_mode == "min") {
    X1_yhat_lasso_train = predict(X1_lasso_fit, train_X1, s = "lambda.min")
    X1_yhat_lasso_test = predict(X1_lasso_fit, test_X1, s = "lambda.min")
    index_X1 = which(X1_lasso_fit$lambda == X1_lasso_fit$lambda.min)
  } else if (fit_mode == "1se") {
    X1_yhat_lasso_train = predict(X1_lasso_fit, train_X1, s = "lambda.1se")
    X1_yhat_lasso_test = predict(X1_lasso_fit, test_X1, s = "lambda.1se")
    index_X1 = which(X1_lasso_fit$lambda == X1_lasso_fit$lambda.1se)
  }
  
  X1_train_e = calc_mse(X1_yhat_lasso_train, train_y)
  X1_err_train_lasso[ii] = X1_train_e
  
  X1_test_e = calc_mse(X1_yhat_lasso_test, test_y)
  X1_err_test_lasso[ii] = X1_test_e
  X1_test_pearson[ii] <- cor(X1_yhat_lasso_test, test_y, method="pearson")
  
  X1_lasso_support[ii] = X1_lasso_fit$nzero[index_X1] 
  print(X1_test_e)
  
  coef_X1_lasso <- coef(X1_lasso_fit, s = "lambda.min", exact = TRUE)
  coef_X1_lasso_dense <- as.matrix(coef_X1_lasso)
  names_X1_lasso <- rownames(coef_X1_lasso_dense)
  coef_X1_lasso_vector <- as.numeric(coef_X1_lasso_dense)
  
  if (names_X1_lasso[1] == "(Intercept)") {
    coef_X1_lasso_vector <- coef_X1_lasso_vector[-1]
    names_X1_lasso <- names_X1_lasso[-1]
  }
  
  # Identify nonzero coefficients
  nonzero_indices <- which(coef_X1_lasso_vector != 0)
  nonzero_X1_lasso <- coef_X1_lasso_vector[nonzero_indices]
  nonzero_names_X1_lasso <- names_X1_lasso[nonzero_indices]
  
  if (length(nonzero_X1_lasso) > 0) {
    temp_data_X1 <- data.frame(iteration=rep(ii, length(nonzero_X1_lasso)),
                               model=rep("X1_lasso", length(nonzero_X1_lasso)),
                               alpha=rep(X1_lasso_fit$lambda.min, length(nonzero_X1_lasso)),
                               coefficient_value=nonzero_X1_lasso,
                               coefficient_name=nonzero_names_X1_lasso)
    nonzero_coeffs <- rbind(nonzero_coeffs, temp_data_X1)
  }
  
  
  # Separate models for X2
  print("Only X2")
  X2_lasso_fit = cv.glmnet(train_X2, train_y,alpha = alpha_value, standardize = F, foldid=foldid)
  
  if (fit_mode == "min") {
    X2_yhat_lasso_train = predict(X2_lasso_fit, train_X2, s = "lambda.min")
    X2_yhat_lasso_test = predict(X2_lasso_fit, test_X2, s = "lambda.min")
    index_X2 = which(X2_lasso_fit$lambda == X2_lasso_fit$lambda.min)
  } else if (fit_mode == "1se") {
    X2_yhat_lasso_train = predict(X2_lasso_fit, train_X2, s = "lambda.1se")
    X2_yhat_lasso_test = predict(X2_lasso_fit, test_X2, s = "lambda.1se")
    index_X2 = which(X2_lasso_fit$lambda == X2_lasso_fit$lambda.1se)
  }
  
  X2_train_e = calc_mse(X2_yhat_lasso_train, train_y)
  X2_err_train_lasso[ii] = X2_train_e
  
  X2_test_e = calc_mse(X2_yhat_lasso_test, test_y)
  X2_err_test_lasso[ii] = X2_test_e
  X2_test_pearson[ii] <- cor(X2_yhat_lasso_test, test_y, method="pearson")
  
  X2_lasso_support[ii] = X2_lasso_fit$nzero[index_X2] 
  print(X2_test_e)
  
  coef_X2_lasso <- coef(X2_lasso_fit, s = "lambda.min", exact = TRUE)
  
  coef_X2_lasso_dense <- as.matrix(coef_X2_lasso)
  names_X2_lasso <- rownames(coef_X2_lasso_dense)
  coef_X2_lasso_vector <- as.numeric(coef_X2_lasso_dense)
  
  if (names_X2_lasso[1] == "(Intercept)") {
    coef_X2_lasso_vector <- coef_X2_lasso_vector[-1]
    names_X2_lasso <- names_X2_lasso[-1]
  }
  
  # Identify nonzero coefficients
  nonzero_indices <- which(coef_X2_lasso_vector != 0)
  nonzero_X2_lasso <- coef_X2_lasso_vector[nonzero_indices]
  nonzero_names_X2_lasso <- names_X2_lasso[nonzero_indices]
  
  if (length(nonzero_X2_lasso) > 0) {
    temp_data_X2 <- data.frame(iteration=rep(ii, length(nonzero_X2_lasso)),
                               model=rep("X2_lasso", length(nonzero_X2_lasso)),
                               alpha=rep(X2_lasso_fit$lambda.min, length(nonzero_X2_lasso)),
                               coefficient_value=nonzero_X2_lasso,
                               coefficient_name=nonzero_names_X2_lasso)
    nonzero_coeffs <- rbind(nonzero_coeffs, temp_data_X2)
  }
  
  # Separate models for X3
  print("Only X3")
  X3_lasso_fit = cv.glmnet(train_X3, train_y,alpha = alpha_value, standardize = F, foldid=foldid)
  
  if (fit_mode == "min") {
    X3_yhat_lasso_train = predict(X3_lasso_fit, train_X3, s = "lambda.min")
    X3_yhat_lasso_test = predict(X3_lasso_fit, test_X3, s = "lambda.min")
    index_X3 = which(X3_lasso_fit$lambda == X3_lasso_fit$lambda.min)
  } else if (fit_mode == "1se") {
    X3_yhat_lasso_train = predict(X3_lasso_fit, train_X3, s = "lambda.1se")
    X3_yhat_lasso_test = predict(X3_lasso_fit, test_X3, s = "lambda.1se")
    index_X3 = which(X3_lasso_fit$lambda == X3_lasso_fit$lambda.1se)
  }
  
  X3_train_e = calc_mse(X3_yhat_lasso_train, train_y)
  X3_err_train_lasso[ii] = X3_train_e
  
  X3_test_e = calc_mse(X3_yhat_lasso_test, test_y)
  X3_err_test_lasso[ii] = X3_test_e
  X3_test_pearson[ii] <- cor(X3_yhat_lasso_test, test_y, method="pearson")
  
  X3_lasso_support[ii] = X3_lasso_fit$nzero[index_X3] 
  print(X3_test_e)
  
  coef_X3_lasso <- coef(X3_lasso_fit, s = "lambda.min", exact = TRUE)
  
  coef_X3_lasso_dense <- as.matrix(coef_X3_lasso)
  names_X3_lasso <- rownames(coef_X3_lasso_dense)
  coef_X3_lasso_vector <- as.numeric(coef_X3_lasso_dense)
  
  if (names_X3_lasso[1] == "(Intercept)") {
    coef_X3_lasso_vector <- coef_X3_lasso_vector[-1]
    names_X3_lasso <- names_X3_lasso[-1]
  }
  
  # Identify nonzero coefficients
  nonzero_indices <- which(coef_X3_lasso_vector != 0)
  nonzero_X3_lasso <- coef_X3_lasso_vector[nonzero_indices]
  nonzero_names_X3_lasso <- names_X3_lasso[nonzero_indices]
  
  if (length(nonzero_X3_lasso) > 0) {
    temp_data_X3 <- data.frame(iteration=rep(ii, length(nonzero_X3_lasso)),
                               model=rep("X3_lasso", length(nonzero_X3_lasso)),
                               alpha=rep(X3_lasso_fit$lambda.min, length(nonzero_X3_lasso)),
                               coefficient_value=nonzero_X3_lasso,
                               coefficient_name=nonzero_names_X3_lasso)
    nonzero_coeffs <- rbind(nonzero_coeffs, temp_data_X3)
  }
  
  
  # Late fusion
  print("Late Fusion")
  second_stage_smp = floor(val_frac * nrow(train_X1))
  val_ind = sort(sample(seq_len(nrow(train_X1)), size = second_stage_smp))
  train_late_ind = setdiff(seq_len(nrow(train_X1)), val_ind)
  
  # Split training and validation sets for each view
  val_X1 = train_X1[val_ind,]
  train_X1_late = train_X1[train_late_ind,]
  
  val_X2 = train_X2[val_ind,]
  train_X2_late = train_X2[train_late_ind,]
  
  val_X3 = train_X3[val_ind,]
  train_X3_late = train_X3[train_late_ind,]
  
  val_y = train_y[val_ind]
  train_y_late = train_y[train_late_ind]
  
  # Fit Lasso models for each view
  X1_lasso_fit_late = cv.glmnet(train_X1_late, train_y_late,alpha = alpha_value, standardize = F)
  X2_lasso_fit_late = cv.glmnet(train_X2_late, train_y_late,alpha = alpha_value, standardize = F)
  X3_lasso_fit_late = cv.glmnet(train_X3_late, train_y_late,alpha = alpha_value, standardize = F)
  
  if (fit_mode == "min") {
    X1_yhat_lasso_late_val = predict(X1_lasso_fit_late, val_X1, s = "lambda.min")
    X1_yhat_lasso_late_test = predict(X1_lasso_fit_late, test_X1, s = "lambda.min")
    X2_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.min")
    X2_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.min")
    X3_yhat_lasso_late_val = predict(X3_lasso_fit_late, val_X3, s = "lambda.min")
    X3_yhat_lasso_late_test = predict(X3_lasso_fit_late, test_X3, s = "lambda.min")
    
    late_index_X1 = which(X1_lasso_fit_late$lambda == X1_lasso_fit_late$lambda.min)
    late_index_X2 = which(X2_lasso_fit_late$lambda == X2_lasso_fit_late$lambda.min)
    late_index_X3 = which(X3_lasso_fit_late$lambda == X3_lasso_fit_late$lambda.min)
  } else if (fit_mode == "1se") {
    X1_yhat_lasso_late_val = predict(X1_lasso_fit_late, val_X1, s = "lambda.1se")
    X1_yhat_lasso_late_test = predict(X1_lasso_fit_late, test_X1, s = "lambda.1se")
    X2_yhat_lasso_late_val = predict(X2_lasso_fit_late, val_X2, s = "lambda.1se")
    X2_yhat_lasso_late_test = predict(X2_lasso_fit_late, test_X2, s = "lambda.1se")
    X3_yhat_lasso_late_val = predict(X3_lasso_fit_late, val_X3, s = "lambda.1se")
    X3_yhat_lasso_late_test = predict(X3_lasso_fit_late, test_X3, s = "lambda.1se")
    
    late_index_X1 = which(X1_lasso_fit_late$lambda == X1_lasso_fit_late$lambda.1se)
    late_index_X2 = which(X2_lasso_fit_late$lambda == X2_lasso_fit_late$lambda.1se)
    late_index_X3 = which(X3_lasso_fit_late$lambda == X3_lasso_fit_late$lambda.1se)
  }
  
  # Create fusion data and fit the fusion model
  fuse_data = data.frame(
    y = val_y,
    X1_pred = as.vector(X1_yhat_lasso_late_val),
    X2_pred = as.vector(X2_yhat_lasso_late_val),
    X3_pred = as.vector(X3_yhat_lasso_late_val)
  )
  
  fit_fuse = lm(y ~ X1_pred + X2_pred + X3_pred, data = fuse_data)
  
  # Predict and calculate MSE for the test set
  fuse_pred_test = predict(fit_fuse, data.frame(
    X1_pred = as.vector(X1_yhat_lasso_late_test),
    X2_pred = as.vector(X2_yhat_lasso_late_test),
    X3_pred = as.vector(X3_yhat_lasso_late_test)
  ))
  
  err_fuse[ii] = calc_mse(fuse_pred_test, test_y)
  print(err_fuse[ii])
  fuse_test_pearson[ii] <- cor(fuse_pred_test, test_y, method="pearson")
  
  
  support_fuse_late[ii] = X1_lasso_fit_late$nzero[late_index_X1] + 
    X2_lasso_fit_late$nzero[late_index_X2] + 
    X3_lasso_fit_late$nzero[late_index_X3]
  
  
  # Early fusion
  print("Early Fusion")
  fit_lasso_cv = cv.glmnet(cbind(train_X1, train_X2, train_X3), train_y,alpha = alpha_value, standardize = F, foldid = foldid)
  
  if (fit_mode == "min") {
    yhat_lasso_train = predict(fit_lasso_cv, cbind(train_X1, train_X2, train_X3), s = "lambda.min")
    yhat_lasso_test = predict(fit_lasso_cv, cbind(test_X1, test_X2, test_X3), s = "lambda.min")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.min)
  } else if (fit_mode == "1se") {
    yhat_lasso_train = predict(fit_lasso_cv, cbind(train_X1, train_X2, train_X3), s = "lambda.1se")
    yhat_lasso_test = predict(fit_lasso_cv, cbind(test_X1, test_X2, test_X3), s = "lambda.1se")
    index_early = which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se)
  }
  
  train_e = calc_mse(yhat_lasso_train, train_y)
  err_train_lasso[ii] = train_e
  lasso_support[ii] = fit_lasso_cv$nzero[index_early]
  test_e = calc_mse(yhat_lasso_test, test_y)
  err_test_lasso[ii] = test_e
  early_fusion_pearson[ii] <- cor(yhat_lasso_test, test_y, method="pearson")
  
  
  
  print(test_e)
  print(lasso_support[ii])
  
  
  #cooperative regression
  print("Cooperative Regression")
  cvm_min = rep(0, times = length(alphalist))
  test_MSE_min = rep(0, times = length(alphalist))
  support_min = rep(0, times = length(alphalist))
  best_fit_coeffs <- list()
  coop_pearson = length(alphalist)
  
  cvm_min_no_pf = rep(0, times = length(alphalist))
  test_MSE_min_no_pf = rep(0, times = length(alphalist))
  support_min_no_pf = rep(0, times = length(alphalist))
  best_fit_coeffs_no_pf <- list()
  coop_no_pf_pearson = length(alphalist)
  
  for(j in 1:length(alphalist)){
    alpha = alphalist[j]
    print(alpha)
    
    pf_null = rep(1, ncol(train_X1) + ncol(train_X2) + ncol(train_X3))
    full_fit_no_pf = coop_cv_multi(train_X1,train_X2,train_X3,train_y,
                                   alpha=alpha,foldid=foldid,
                                   nfolds=max(foldid),
                                   pf_values=pf_null,
                                   fit_mode = fit_mode)
    
    yhat_coop_new_no_pf = cbind(train_X1,train_X2,train_X3) %*% full_fit_no_pf$best_fit_coef + (full_fit_no_pf$best_fit_intercept)
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, train_y)
    
    if (fit_mode == 'min'){
      new_cv_coop_no_pf[ii, j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_min]
      cvm_min_no_pf[j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_min]
      new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[full_fit_no_pf$ind_min]
    } else if (fit_mode == '1se'){
      new_cv_coop_no_pf[ii, j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_1se]
      cvm_min_no_pf[j] = full_fit_no_pf$cvm[full_fit_no_pf$ind_1se]
      new_support_coop_no_pf[ii, j] = full_fit_no_pf$support[full_fit_no_pf$ind_1se]
    }
    
    new_train_coop_no_pf[ii,j] = calc_mse(yhat_coop_new_no_pf, train_y)
    yhat_coop_new_test_no_pf = cbind(test_X1, test_X2, test_X3) %*% full_fit_no_pf$best_fit_coef + (full_fit_no_pf$best_fit_intercept)
    test_e_coop_new_no_pf = calc_mse(yhat_coop_new_test_no_pf, test_y)
    new_test_coop_no_pf[ii,j] = test_e_coop_new_no_pf
    coop_no_pf_pearson[j] <- cor(as.matrix(yhat_coop_new_test_no_pf), test_y, method="pearson")
    
    print("Full (no penalty factor)")
    print(new_test_coop_no_pf[ii,j])
    
    test_MSE_min_no_pf[j] = test_e_coop_new_no_pf
    support_min_no_pf[j] = new_support_coop_no_pf[ii, j]
    
    best_fit_coeffs_no_pf[[paste("Alpha", alpha)]] <- as.matrix(full_fit_no_pf$best_fit_coef)
    
    
    #Iterative Coop Learning
    coop_fit_iter = coop_regression_iter_multi(train_X1,train_X2,train_X3,train_y,
                                               alpha=alpha,foldid=foldid,
                                               max_iteration=5,
                                               thr=0.05)
    yhat_coop_iter = (train_X1%*%coop_fit_iter$thetax1 + coop_fit_iter$intercept_x1) +
      (train_X2%*%coop_fit_iter$thetax2 + coop_fit_iter$intercept_x2) +
      (train_X3%*%coop_fit_iter$thetax3 + coop_fit_iter$intercept_x3)
    support_coop[ii, j] = coop_fit_iter$n_thetax1 + coop_fit_iter$n_thetax2 + coop_fit_iter$n_thetax3
    
    yhat_coop_test = (test_X1%*%coop_fit_iter$thetax1 + coop_fit_iter$intercept_x1) +
      (test_X2%*%coop_fit_iter$thetax2 + coop_fit_iter$intercept_x2) +
      (test_X3%*%coop_fit_iter$thetax3 + coop_fit_iter$intercept_x3)
    test_e_coop = calc_mse(yhat_coop_test, test_y)
    err_test_coop[ii,j] = test_e_coop
    print("Iterative")
    
    #Iterative Coop Learning
    lambda_x1 = coop_fit_iter$lam_x1
    lambda_x2 = coop_fit_iter$lam_x2
    lambda_x3 = coop_fit_iter$lam_x3
    
    nx1 = ncol(train_X1)
    nx2 = ncol(train_X2)
    nx3 = ncol(train_X3)
    
    adjust_factor = (nx1 * lambda_x1 + nx2 * lambda_x2 + + nx3 * lambda_x3) / (nx1 + nx2 + nx3)
    pfs = c(rep(lambda_x1, ncol(train_X1)),rep(lambda_x2, ncol(train_X2)),rep(lambda_x3, ncol(train_X3)))
    full_fit = coop_cv_multi(train_X1,train_X2,train_X3,train_y,
                             alpha=alpha,foldid=foldid,
                             nfolds=max(foldid),
                             pf_values=pfs,
                             fit_mode = fit_mode)
    
    yhat_coop_new = cbind(train_X1,train_X2,train_X3) %*% full_fit$best_fit_coef + (full_fit$best_fit_intercept)
    new_train_coop[ii,j] = calc_mse(yhat_coop_new, train_y)
    
    if (fit_mode == 'min'){
      new_cv_coop[ii, j] = full_fit$cvm[full_fit$ind_min]
      cvm_min[j] = full_fit$cvm[full_fit$ind_min]
      new_support_coop[ii, j] = full_fit$support[full_fit$ind_min]
    } else if (fit_mode == '1se'){
      new_cv_coop[ii, j] = full_fit$cvm[full_fit$ind_1se]
      cvm_min[j] = full_fit$cvm[full_fit$ind_1se]
      new_support_coop[ii, j] = full_fit$support[full_fit$ind_1se]
    }
    
    new_train_coop[ii,j] = calc_mse(yhat_coop_new, train_y)
    yhat_coop_new_test = cbind(test_X1, test_X2, test_X3) %*% full_fit$best_fit_coef + (full_fit$best_fit_intercept)
    test_e_coop_new = calc_mse(yhat_coop_new_test, test_y)
    new_test_coop[ii,j] = test_e_coop_new
    coop_pearson[j] <- cor(as.matrix(yhat_coop_new_test), test_y, method="pearson")
    
    print("Full (with penalty factor)")
    print(new_test_coop[ii,j])
    
    test_MSE_min[j] = test_e_coop_new
    support_min[j] = new_support_coop[ii, j]
    best_fit_coeffs[[paste("Alpha", alpha)]] <- as.matrix(full_fit$best_fit_coef)
    
  }
  
  coop_selected_by_cv_no_pf[ii] = test_MSE_min_no_pf[which.min(cvm_min_no_pf)]
  support_by_cv_no_pf[ii] = support_min_no_pf[which.min(cvm_min_no_pf)]
  alpha_by_cv_no_pf[ii] = alphalist[which.min(cvm_min_no_pf)]
  coop_coeffs_no_pf[ii] = best_fit_coeffs_no_pf[which.min(cvm_min_no_pf)]
  coop_no_pf_pearson_best[ii] = coop_no_pf_pearson[which.min(cvm_min_no_pf)]
  
  
  print("Full selected by cv, without pf")
  print(coop_selected_by_cv_no_pf[ii])
  print(alpha_by_cv_no_pf[ii])
  
  coop_selected_by_cv[ii] = test_MSE_min[which.min(cvm_min)]
  support_by_cv[ii] = support_min[which.min(cvm_min)]
  alpha_by_cv[ii] = alphalist[which.min(cvm_min)]
  coop_coeffs[ii] = best_fit_coeffs[which.min(cvm_min)]
  coop_pearson_best[ii] = coop_pearson[which.min(cvm_min)]
  
  
  
  print("Full selected by cv, with pf")
  print(coop_selected_by_cv[ii])
  print(alpha_by_cv[ii])
}


test_err = sqrt(cbind(err_null,X1_err_test_lasso, X2_err_test_lasso, X3_err_test_lasso, 
                 err_test_lasso, err_fuse, coop_selected_by_cv_no_pf, 
                 coop_selected_by_cv, new_test_coop_no_pf, new_test_coop))

support_df = cbind(X1_lasso_support, X2_lasso_support, X3_lasso_support, 
                   support_fuse_late, lasso_support, support_by_cv_no_pf, 
                   support_by_cv, new_support_coop_no_pf, new_support_coop)

pearson_df = cbind(X1_test_pearson, X2_test_pearson, X3_test_pearson,
                   fuse_test_pearson, early_fusion_pearson,
                   coop_no_pf_pearson_best, coop_pearson_best)

out = rbind(colMeans(test_err), apply(test_err, 2, sd) / sqrt(simN))
test_err_all_diff = test_err - test_err[, 4]  
out_diff = rbind(colMeans(test_err_all_diff), apply(test_err_all_diff, 2, sd) / sqrt(10))
colMeans(support_df[, 1:9])  


saveRDS(test_err, paste0("meal_test_err.rds"))
saveRDS(support_df, paste0("meal_support.rds"))
saveRDS(alpha_by_cv_no_pf, paste0("meal_alpha_no_pf.rds"))
saveRDS(alpha_by_cv, paste0("meal_alpha_pf.rds"))
write.csv(nonzero_coeffs, "nonzero_coefficients.csv", row.names=FALSE)
saveRDS(coop_coeffs, paste0("coop_coeffs.rds"))
saveRDS(coop_coeffs_no_pf, paste0("coop_coeffs_no_pf.rds"))
saveRDS(pearson_df, paste0("pearson_df.rds"))


