library(ggplot2); library(sandwich)

make_lin_DGP = function(n, t, ATE, beta_x, seed, sd = 1) {
  set.seed(seed)
  trt = sample(rep(c(0, 1), c(n/2, n/2)), size = n, replace = FALSE)
  
  X = rnorm(n, mean = 1, sd = sd)
  X[trt == 1] = X[trt == 0]
  
  y = list()
  for (i in 1:t) {
    y[[i]] = rnorm(n) + beta_x*X 
    y[[i]][trt == 1] = rnorm(sum(trt == 1), mean = ATE) + beta_x*X[trt == 1]
    
  }
  
  l = list(trt, y, X)
  return(l)
}

make_nonlin_DGP = function(n, t, ATE, beta_x, seed, sd = 0.5) {
  set.seed(seed)
  trt = sample(rep(c(0, 1), c(n/2, n/2)), size = n, replace = FALSE)
  
  X = rnorm(n, mean = 1, sd = sd)
  X[trt == 1] = X[trt == 0]
  
  y = list()
  for (i in 1:t) {
    y[[i]] = rexp(n) + beta_x*X^2*log(abs(X)) 
    y[[i]][trt == 1] = rexp(sum(trt == 1)) + beta_x*X[trt == 1]^2*log(abs(X[trt == 1]))  + ATE
    
  }
  
  l = list(trt, y, X)
  return(l)
}

make_nonlin_DGP2 = function(n, t, ATE, beta_x, seed, sd = 0.5) {
  set.seed(seed)
  trt = sample(rep(c(0, 1), c(n/2, n/2)), size = n, replace = FALSE)
  
  X = rnorm(n, mean = 1, sd = sd)
  X[trt == 1] = X[trt == 0]
  
  y = list()
  for (i in 1:t) {
    y[[i]] = rexp(n) + beta_x*X^2*log(abs(X)) 
    y[[i]][trt == 1] = rexp(sum(trt == 1)) + beta_x*X[trt == 1]^2*log(abs(X[trt == 1])) + ATE/(1 + log(i))
    
  }
  
  l = list(trt, y, X)
  return(l)
}

get_CS_no_cov = function(l) {
  y = l[[2]]; trt = l[[1]]; n = length(y[[1]]); t = length(y)
  
  trt_df = ctrl_df = data.frame(NA)
  hatbeta = hatsigma = vector()
  for (i in 1:t) {
    in_y = y[[i]]
    in_y_trt = y[[i]][trt == 1]
    in_y_ctrl = y[[i]][trt == 0]
    trt_df = cbind(trt_df, in_y_trt)
    ctrl_df = cbind(ctrl_df, in_y_ctrl)
    lm_mod = lm(in_y ~ trt)
    hatbeta[i] = coef(lm_mod)[2]
    hatsigma[i] = as.numeric(diag(vcovHC(lm_mod, type = "HC")))[2]
  }
  
  ## acounting for lagged dependency
  trt_df = trt_df[, -1]; ctrl_df = ctrl_df[, -1]
  
  Shat1 = cov(trt_df)
  Shat0 = cov(ctrl_df)
  
  n1 = sum(trt == 1); n0 = sum(trt == 0)
  
  sigma_t1 = Shat1[1, 1]
  sigma_t0 = Shat0[1, 1]
  s = sigma_t1/n1 + sigma_t0/n0
  for (i in 2:t) {
    sigma_t1[i] = (Shat1[i, i] - as.numeric(matrix(Shat1[i, (1:(i-1))], ncol = (i- 1))%*% solve(Shat1[1:(i-1), 1:(i-1)]) %*% matrix(Shat1[(1:(i-1)), i], nrow = (i-1))))/n1
    sigma_t0[i] = (Shat0[i, i] - as.numeric(matrix(Shat0[i, (1:(i-1))], ncol = (i- 1))%*% solve(Shat0[1:(i-1), 1:(i-1)]) %*% matrix(Shat0[(1:(i-1)), i], nrow = (i-1))))/n0
    s[i] = sigma_t1[i] + sigma_t0[i]
  }
  
  df = data.frame(center = hatbeta, adjusted_var = s, unadjusted_var = hatsigma)
  return(df)
}

get_CS_w_cov = function(l) {
  y = l[[2]]; trt = l[[1]]; n = length(y[[1]]); t = length(y); X = l[[3]]
  
  trt_df = ctrl_df = data.frame(NA)
  hatbeta = hatsigma = vector()
  hatbeta_proxy = vector()
  for (i in 1:t) {
    in_y = y[[i]]
    lm_mod = lm(in_y ~ trt + X)
    hatbeta[i] = coef(lm_mod)[2]
    hatsigma[i] = as.numeric(diag(vcovHC(lm_mod, type = "HC")))[2]
    
    resid_lm_mod = lm(in_y ~ X)
    new_y = residuals(resid_lm_mod)
    in_y_trt = new_y[trt == 1]
    in_y_ctrl = new_y[trt == 0]
    trt_df = cbind(trt_df, in_y_trt)
    ctrl_df = cbind(ctrl_df, in_y_ctrl)
    hatbeta_proxy[i] = mean(in_y_trt) - mean(in_y_ctrl)
    
  }
  
  ## acounting for lagged dependency
  
  
  
  trt_df = trt_df[, -1]; ctrl_df = ctrl_df[, -1]
  Shat1 = cov(trt_df)
  Shat0 = cov(ctrl_df)
  
  n1 = sum(trt == 1); n0 = sum(trt == 0)
  
  sigma_t1 = Shat1[1, 1]
  sigma_t0 = Shat0[1, 1]
  s = sigma_t1/n1 + sigma_t0/n0
  for (i in 2:t) {
    sigma_t1[i] = (Shat1[i, i] - as.numeric(matrix(Shat1[i, (1:(i-1))], ncol = (i- 1))%*% solve(Shat1[1:(i-1), 1:(i-1)]) %*% matrix(Shat1[(1:(i-1)), i], nrow = (i-1))))/n1
    sigma_t0[i] = (Shat0[i, i] - as.numeric(matrix(Shat0[i, (1:(i-1))], ncol = (i- 1))%*% solve(Shat0[1:(i-1), 1:(i-1)]) %*% matrix(Shat0[(1:(i-1)), i], nrow = (i-1))))/n0
    s[i] = sigma_t1[i] + sigma_t0[i]
  }
  
  df = data.frame(center = hatbeta, unadjusted_var = hatsigma, adjusted_center = hatbeta_proxy, adjusted_var = s)
  return(df)
}

get_CS = function(hatbeta, s, eta = 0.77, alpha = 0.05) {
  t = length(s)
  S = cumsum(s)
  center = cumsum(hatbeta)/(1:t)
  radius = sqrt((S*eta^2 + 1)/eta^2 * log((S*eta^2 + 1)/alpha^2))*(1/(1:t))
  lower = center - radius
  upper = center + radius
  df = data.frame(center, lower, upper)
  return(df)
}

get_stats = function(get_CS_obj, truth, check = 1, in_t) {
  stoppingtime = which(get_CS_obj[, 2] >= 0)[1]
  width = get_CS_obj[in_t, 3] - get_CS_obj[in_t, 2]
  all_checks = ((get_CS_obj[, 2] <= truth) & (get_CS_obj[, 3] >= truth))
  subset_checks = all_checks[check:length(all_checks)]
  
  coverage = all(subset_checks)
  result_df = data.frame(stoppingtime, width, coverage)
  return(result_df)
}

n = 1000; t = 15; ATE = 0.5; beta_x = 5.0; check = 1; eta = 3.0; in_t = 15; sd = 1
no_cov_unadj_df = no_cov_adj_df = data.frame()
cov_unadj_df = cov_adj_df = data.frame()

B = 1000
for (i in 1:B) {
  l = make_lin_DGP(n = n, t = t, ATE = ATE, beta_x = beta_x, seed = i, sd = sd)
  
  no_cov = get_CS_no_cov(l)
  
  no_cov_unadj = get_CS(no_cov[, 1], s = no_cov$unadjusted_var, eta = eta)
  no_cov_unadj_result = get_stats(no_cov_unadj, truth = ATE, check = check, in_t = in_t)
  no_cov_unadj_df = rbind(no_cov_unadj_df, no_cov_unadj_result)
  
  no_cov_adj = get_CS(no_cov[, 1], s = no_cov$adjusted_var, eta = eta)
  no_cov_adj_result = get_stats(no_cov_adj, truth = ATE, check = check, in_t = in_t)
  no_cov_adj_df = rbind(no_cov_adj_df, no_cov_adj_result)
  
  cov = get_CS_w_cov(l)
  
  cov_unadj = get_CS(cov[, 1], s = cov$unadjusted_var, eta = eta)
  cov_unadj_result = get_stats(cov_unadj, truth = ATE, check = check, in_t = in_t)
  cov_unadj_df = rbind(cov_unadj_df, cov_unadj_result)
  
  cov_adj = get_CS(cov[, 3], s = cov$adjusted_var, eta = eta)
  cov_adj_result = get_stats(cov_adj, truth = ATE, check = check, in_t = in_t)
  cov_adj_df = rbind(cov_adj_df, cov_adj_result)
  
  print(i)
}

apply(no_cov_unadj_df, 2, mean); apply(no_cov_adj_df, 2, mean)
apply(cov_unadj_df, 2, mean); apply(cov_adj_df, 2, mean)


####
n = 1000; t = 30; ATE = 0.5; beta_x = 5.0; check = 1; eta = 3.0; in_t = 15; sd = 1.0; t2 = 70
no_cov_unadj_df = no_cov_adj_df = data.frame()
cov_unadj_df = cov_adj_df = data.frame()

B = 1000
for (i in 1:B) {
  l = make_nonlin_DGP(n = n, t = t2, ATE = ATE, beta_x = beta_x, seed = i, sd = sd)
  
  no_cov = get_CS_no_cov(l)
  
  no_cov_unadj = get_CS(no_cov[, 1], s = no_cov$unadjusted_var, eta = eta)
  no_cov_unadj_result = get_stats(no_cov_unadj, truth = ATE, check = check, in_t = in_t)
  no_cov_unadj_df = rbind(no_cov_unadj_df, no_cov_unadj_result)
  
  no_cov_adj = get_CS(no_cov[, 1], s = no_cov$adjusted_var, eta = eta)
  no_cov_adj_result = get_stats(no_cov_adj, truth = ATE, check = check, in_t = in_t)
  no_cov_adj_df = rbind(no_cov_adj_df, no_cov_adj_result)
  
  cov = get_CS_w_cov(l)
  
  cov_unadj = get_CS(cov[, 1], s = cov$unadjusted_var, eta = eta)
  cov_unadj_result = get_stats(cov_unadj, truth = ATE, check = check, in_t = in_t)
  cov_unadj_df = rbind(cov_unadj_df, cov_unadj_result)
  
  cov_adj = get_CS(cov[, 3], s = cov$adjusted_var, eta = eta)
  cov_adj_result = get_stats(cov_adj, truth = ATE, check = check, in_t = in_t)
  cov_adj_df = rbind(cov_adj_df, cov_adj_result)
  print(i)
}

apply(no_cov_unadj_df, 2, mean); apply(no_cov_adj_df, 2, mean)
apply(cov_unadj_df, 2, mean); apply(cov_adj_df, 2, mean)
summary(no_cov_unadj_df$stoppingtime) #36.9

######
n = 1000; t = 30; ATE = 2; beta_x = 5.0; check = 1; eta = 3.0; in_t = 15; sd = 1.0; t2 = 40
no_cov_unadj_df = no_cov_adj_df = data.frame()
cov_unadj_df = cov_adj_df = data.frame()
truth = ATE/(1 + log(1:t2))
truth = cumsum(truth)/(1:t2)
B = 50
for (i in 1:B) {
  l = make_nonlin_DGP2(n = n, t = t2, ATE = ATE, beta_x = beta_x, seed = i, sd = sd)
  
  no_cov = get_CS_no_cov(l)
  
  no_cov_unadj = get_CS(no_cov[, 1], s = no_cov$unadjusted_var, eta = eta)
  no_cov_unadj_result = get_stats(no_cov_unadj, truth = truth, check = check, in_t = in_t)
  no_cov_unadj_df = rbind(no_cov_unadj_df, no_cov_unadj_result)
  
  no_cov_adj = get_CS(no_cov[, 1], s = no_cov$adjusted_var, eta = eta)
  no_cov_adj_result = get_stats(no_cov_adj, truth = truth, check = check, in_t = in_t)
  no_cov_adj_df = rbind(no_cov_adj_df, no_cov_adj_result)
  
  cov = get_CS_w_cov(l)
  
  cov_unadj = get_CS(cov[, 1], s = cov$unadjusted_var, eta = eta)
  cov_unadj_result = get_stats(cov_unadj, truth = truth, check = check, in_t = in_t)
  cov_unadj_df = rbind(cov_unadj_df, cov_unadj_result)
  
  cov_adj = get_CS(cov[, 3], s = cov$adjusted_var, eta = eta)
  cov_adj_result = get_stats(cov_adj, truth = truth, check = check, in_t = in_t)
  cov_adj_df = rbind(cov_adj_df, cov_adj_result)
  print(i)
}

apply(no_cov_unadj_df, 2, mean); apply(no_cov_adj_df, 2, mean)
apply(cov_unadj_df, 2, mean); apply(cov_adj_df, 2, mean)
summary(no_cov_unadj_df$stoppingtime) #36.9

### just for lower n asymptotics taking longer
n = 100; t = 30; ATE = 0.5; beta_x = 0; check = 1; eta = 3.0; in_t = 15; sd = 1.0; t2 = t
no_cov_unadj_df = no_cov_adj_df = data.frame()
cov_unadj_df = cov_adj_df = data.frame()

B = 1000
for (i in 1:B) {
  l = make_nonlin_DGP(n = n, t = t2, ATE = ATE, beta_x = beta_x, seed = i, sd = sd)
  
  no_cov = get_CS_no_cov(l)
  
  no_cov_unadj = get_CS(no_cov[, 1], s = no_cov$unadjusted_var, eta = eta)
  no_cov_unadj_result = get_stats(no_cov_unadj, truth = ATE, check = check, in_t = in_t)
  no_cov_unadj_df = rbind(no_cov_unadj_df, no_cov_unadj_result)
  
  no_cov_adj = get_CS(no_cov[, 1], s = no_cov$adjusted_var, eta = eta)
  no_cov_adj_result = get_stats(no_cov_adj, truth = ATE, check = check, in_t = in_t)
  no_cov_adj_df = rbind(no_cov_adj_df, no_cov_adj_result)
  print(i)
}

apply(no_cov_unadj_df, 2, mean); apply(no_cov_adj_df, 2, mean)




