.libPaths(c("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackage_rstan/"))
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(rstan)
library(truncnorm)
library(stringr)
set.seed(123)


case_id_list <- c("ABNCR6D")
  #c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

# ------------------------------------------------------------------------------
# MCMC model
for(case_id in case_id_list){
  print(case_id)
  
  cell_cover_germline_mutation_site <- readRDS(paste0("./results/germline_mutation_VAF/", case_id, ".rds"))
  cell_cover_germline_mutation_site <- na.omit(cell_cover_germline_mutation_site)
  cell_cover_germline_mutation_site <- cell_cover_germline_mutation_site[cell_cover_germline_mutation_site$covered > 0,]
  # summary(lm(detected ~ covered, cell_cover_germline_mutation_site))
  
  # Define the Stan model
  stan_code <- "
  data {
    int<lower=0> N;         // Number of data points
    vector<lower=0>[N] x;            // Predictor variable
    vector<lower=0>[N] y;            // Response variable
  }
  parameters {
    real alpha;             // Intercept
    real<lower=0, upper=1> beta;              // Slope
    real<lower=0> sigma;    // Standard deviation of the errors
  }
  model {
    // Prior distributions
    alpha ~ normal(0, 100);
    beta ~ normal(0.5, 10000);
    sigma ~ cauchy(0, 10000);
    
    // Likelihood
      y ~ normal(alpha + beta * x, sigma);
  }
"
  
  # Compile the Stan model
  stan_model <- stan_model(model_code = stan_code)
  # Prepare data list for Stan
  stan_data <- list(N = dim(cell_cover_germline_mutation_site)[1], 
                    x = cell_cover_germline_mutation_site$covered, 
                    y = cell_cover_germline_mutation_site$detected)
  # Fit the Bayesian linear regression model using MCMC
  fit <- sampling(stan_model, data = stan_data, chains = 4, iter = 2000)
  print(fit)
  saveRDS(fit, paste0("./results/germline_mutation_VAF/MCMC_gaussian/", case_id, ".rds"))
}


# ------------------------------------------------------------------------------
# posterior distribution
for(case_id in case_id_list){
  print(case_id)
  
  # Read in clonal mutations
  cell_cover_clonal_mutation_site <- 
    readRDS(paste0("./results/clonal_mutation/clonal_mutation_VAF_", case_id, ".rds"))
  
  if(dim(cell_cover_clonal_mutation_site)[1] == 0){next}
  
  # Extract samples from the MCMC output
  fit <- readRDS(paste0("./results/germline_mutation_VAF/MCMC_gaussian/", case_id, ".rds"))
  trace <- extract(fit)
  
  # Extract parameters from trace
  alpha_trace <- trace$alpha
  beta_trace <- trace$beta
  sigma_trace <- trace$sigma
  
  # Calculate posterior distribution 
  for (i in 1:dim(cell_cover_clonal_mutation_site)[1]){
    x = cell_cover_clonal_mutation_site[i, "covered"]
    posterior_y <- rtruncnorm(nrow(alpha_trace), a = 0, mean = alpha_trace + beta_trace * x, sd = sigma_trace)
    
    # Visualize the posterior distribution of y
    # hist(posterior_y, main = "Posterior Distribution of y for Given x", xlab = "y")
    
    cell_cover_clonal_mutation_site[i, "bound_0.001%"] <- quantile(posterior_y, probs = 0.00001)
    cell_cover_clonal_mutation_site[i, "bound_0.01%"] <- quantile(posterior_y, probs = 0.0001)
    cell_cover_clonal_mutation_site[i, "bound_0.1%"] <- quantile(posterior_y, probs = 0.001)
    cell_cover_clonal_mutation_site[i, "bound_1%"] <- quantile(posterior_y, probs = 0.01)
    cell_cover_clonal_mutation_site[i, "bound_5%"] <- quantile(posterior_y, probs = 0.05)
    cell_cover_clonal_mutation_site[i, "bound_95%"] <- quantile(posterior_y, probs = 0.95)
    cell_cover_clonal_mutation_site[i, "bound_99%"] <- quantile(posterior_y, probs = 0.99)
  }
  
  saveRDS(cell_cover_clonal_mutation_site, paste0("./results/clonal_mutation/clonal_mutation_VAF_", case_id, ".rds"))
}

