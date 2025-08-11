library(mlogit)
library(dplyr)
library(readr)     

## Network size
n <- 100

## Sparsity indices
a <- 0.5
b <- 0.5

## True parameters
mu_dagger <- 0.5
rho_dagger <- 0.5
gamma_1 <- 0.2
gamma_2 <- 0.4
delta <- 0.3

mu  <- -a * log(n) + mu_dagger
rho <-  b * log(n) + rho_dagger
theta_0 <- c(mu, rho, gamma_1, gamma_2, delta)

## Convert two dimensional index to one dimension
edge_index <- function(i, j, n) {
  (j - i) + (i - 1) * (2 * n - i) / 2
}

## Number of iterations
M <- 1000

for (rep in 1:M) {
  
  # Nodewise covariates and edgewise covariates
  X <- runif(n, -1, 1)
  Y <- runif(n, -1, 1)
  V <- runif(n * (n - 1) / 2, -1, 1)  
  
  m <- n * (n - 1) / 2
  node_cov <- 2
  edge_cov <- 1
  data_mat <- matrix(0, nrow = m, ncol = 2*(1+node_cov)+edge_cov) 
  
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      idx <- edge_index(i, j, n)
      Vij <- V[idx]
      
      k0 <- 1 + 
        exp(mu + X[i] * gamma_1 + Y[j] * gamma_2) +
        exp(mu + X[j] * gamma_1 + Y[i] * gamma_2) +
        exp(2 * mu + rho + Vij * delta + (X[i] + X[j]) * gamma_1 + (Y[i] + Y[j]) * gamma_2)
      
      probs <- c(
        p00 = 1 / k0,
        p10 = exp(mu + X[i] * gamma_1 + Y[j] * gamma_2) / k0,
        p01 = exp(mu + X[j] * gamma_1 + Y[i] * gamma_2) / k0,
        p11 = exp(2 * mu + rho + Vij * delta + (X[i] + X[j]) * gamma_1 + (Y[i] + Y[j]) * gamma_2) / k0
      )
      
      A <- sample(c("00", "10", "01", "11"), size = 1, prob = probs)
      edge_data <- switch(A,
                          "00" = c(0, 0, X[i], X[j], Y[i], Y[j], Vij),
                          "10" = c(1, 0, X[i], X[j], Y[i], Y[j], Vij),
                          "01" = c(0, 1, X[i], X[j], Y[i], Y[j], Vij),
                          "11" = c(1, 1, X[i], X[j], Y[i], Y[j], Vij)
      )
      
      data_mat[idx, ] <- edge_data
    }
  }
  df <- as.data.frame(data_mat) %>%
    rename(
      A_ij = V1, A_ji = V2,
      X_i = V3, X_j = V4,
      Y_i = V5, Y_j = V6,
      V_edge = V7
    )
  
  df <- df %>%
    mutate(
      choice = case_when(
        A_ij == 0 & A_ji == 0 ~ 1,
        A_ij == 1 & A_ji == 0 ~ 2,
        A_ij == 0 & A_ji == 1 ~ 3,
        TRUE ~ 4
      ),
      mu1 = 0, mu2 = 1, mu3 = 1, mu4 = 2,
      rho1 = 0, rho2 = 0, rho3 = 0, rho4 = 1,
      X1 = 0, X2 = X_i, X3 = X_j, X4 = X_i + X_j,
      Y1 = 0, Y2 = Y_j, Y3 = Y_i, Y4 = Y_i + Y_j,
      Z1 = 0, Z2 = 0, Z3 = 0, Z4 = V_edge
    )
  
  n_dim <- 2 + node_cov + edge_cov
  end <- dim(df)[2]
  start <- dim(df)[2] - n_dim*4 + 1
  mlogit_df <- mlogit.data(
    df,
    choice = "choice",
    shape = "wide",
    varying = start:end,
    sep = ""
  )
  
  model <- mlogit(choice ~ mu + rho + X + Y + Z | 0, data = mlogit_df)
  model_summary <- summary(model)
  for (j in 1:n_dim) {
    estimator <- c(model_summary$CoefTable[j,1],
                   model_summary$CoefTable[j,1]-1.96*model_summary$CoefTable[j,2],
                   model_summary$CoefTable[j,1]+1.96*model_summary$CoefTable[j,2],
                   (theta_0[j] > model_summary$CoefTable[j,1]-1.96*model_summary$CoefTable[j,2])&(theta_0[j] <model_summary$CoefTable[j,1]+1.96*model_summary$CoefTable[j,2]),
                   2*1.96*model_summary$CoefTable[j,2])
    file_name <- paste("output",as.character(j),"_",as.character(n),"_",as.character(a),"_",as.character(b),'.csv', sep="")
    write.table(t(estimator),file = file_name,append = TRUE,
                sep = ",",                 
                row.names = FALSE,          
                col.names = !file.exists(file_name),
                quote = FALSE)
  }
}

par(mfrow = c(1, n_dim),             
    cex.axis = 1.4,               
    cex.lab = 1.4,                
    font.axis = 1.2,              
    font.lab = 1.2,               
    lwd = 2,                      
    mar = c(4, 4, 2, 1))          

for (i in 1:n_dim) {
  file_name <- paste("output",as.character(i),"_",as.character(n),"_",as.character(a),"_",as.character(b),'.csv', sep="")
  data <- read_csv(files)
  stat <- (data$V1 - theta_0[i]) / (data$V5 / (2 * 1.96))
  qqnorm(stat,
         pch = 1,
         cex = 0.85,
         frame = FALSE,
         main = '',
         xlim = c(-3, 3),
         ylim = c(-3, 3),
         xlab = "",
         ylab = "")
  qqline(stat, col = "steelblue", lwd = 2)
}
