library(mlogit)
library(matrixStats)
library(tidyverse)
library(igraph)

library(readxl)

###### Analysis of lawyer friendship network data

attributes <- read.delim("LazegaLawyers/ELattr.dat", 
                         header = FALSE, sep = " ")
names(attributes) <- c("seniority",
                       "status",
                       "gender",
                       "office",
                       "years_with_firm",
                       "age",
                       "practice",
                       "law_school")

# remove unnecessary columns
attributes <- attributes[, -c(1,9)]

adj_adv <- read.delim("LazegaLawyers/ELadv.dat",
                      header = FALSE, sep = " ")
adj_adv <- as.matrix(adj_adv)
adj_adv <- adj_adv[,-1]

## Not symmetric
isSymmetric(adj_adv) 

### summary statistics
max(rowSums(adj_adv)) #30
min(rowSums(adj_adv))  #0

max(colSums(adj_adv)) #37
min(colSums(adj_adv)) #0

mean(colSums(adj_adv)) #mean degree 12.56
mean(colSums(adj_adv))/70 #density 0.18


node <- dim(adj_adv)[1]
dim_X <- 1
dim_Y <- 1
dim_Z <- 5

n_col <- 2 + 2*(dim_X+dim_Y)+ dim_Z
n_row <- node*(node-1)/2

data <- matrix(rep(0,n_row*n_col),nrow = n_row)

for (i in c(1:(node-1))) {
  for (j in c((i+1):node)){
    A_ij <- adj_adv[i,j]
    A_ji <- adj_adv[j,i]
    X_i <- attributes[i,5]
    X_j <- attributes[j,5]
    Y_i <- attributes[i,4]
    Y_j <- attributes[j,4]
    Z_ij <- c(attributes[i,1]==attributes[j,1],
              attributes[i,2]==attributes[j,2],
              attributes[i,3]==attributes[j,3],
              attributes[i,6]==attributes[j,6],
              attributes[i,7]==attributes[j,7])
    data[(j-i)+(i-1)*(2*node-i)/2,] <- c(A_ij,A_ji,X_i,X_j,Y_i,Y_j,Z_ij)
  }
}


df1 <- as.data.frame(data) %>%
  rename(
    A_ij = V1, A_ji = V2,
    X_i = V3, X_j = V4,
    Y_i = V5, Y_j = V6,
    Z1 = V7, Z2 = V8, Z3 = V9, Z4 = V10, Z5 = V11
  )


X_name <- rep(0,4)
Y_name <- rep(0,4)
Z_name <- list()
for (i in 1:4){
  X_name[i] <- paste0(names(attributes)[5],i)
  Y_name[i] <- paste0(names(attributes)[4],i)
}
for (j in 1:3){
  Z <- rep(0,4)
  for (i in 1:4){
    Z[i] <- paste0(names(attributes)[j],i)
  }
  Z_name <- append(Z_name, list(Z))
}

for (j in 1:2){
  Z <- rep(0,4)
  for (i in 1:4){
    Z[i] <- paste0(names(attributes)[j+5],i)
  }
  Z_name <- append(Z_name, list(Z))
}

df1 <- df1 %>%
  mutate(
    choice = case_when(
      A_ij == 0 & A_ji == 0 ~ 1,
      A_ij == 1 & A_ji == 0 ~ 2,
      A_ij == 0 & A_ji == 1 ~ 3,
      TRUE ~ 4
    ),
    mu1 = 0, mu2 = 1, mu3 = 1, mu4 = 2,
    rho1 = 0, rho2 = 0, rho3 = 0, rho4 = 1,
    !!X_name[1] := 0,!!X_name[2] := X_i, !!X_name[3] := X_j, !!X_name[4] := X_i + X_j,
    !!Y_name[1] := 0, !!Y_name[2] := Y_j, !!Y_name[3] := Y_i, !!Y_name[4] := Y_i + Y_j,
    !!Z_name[1][[1]][1] := 0, !!Z_name[1][[1]][2] := 0, !!Z_name[1][[1]][3] := 0, !!Z_name[1][[1]][4] := Z1,
    !!Z_name[2][[1]][1] := 0, !!Z_name[2][[1]][2] := 0, !!Z_name[2][[1]][3] := 0, !!Z_name[2][[1]][4] := Z2,
    !!Z_name[3][[1]][1] := 0, !!Z_name[3][[1]][2] := 0, !!Z_name[3][[1]][3] := 0, !!Z_name[3][[1]][4] := Z3,
    !!Z_name[4][[1]][1] := 0, !!Z_name[4][[1]][2] := 0, !!Z_name[4][[1]][3] := 0, !!Z_name[4][[1]][4] := Z4,
    !!Z_name[5][[1]][1] := 0, !!Z_name[5][[1]][2] := 0, !!Z_name[5][[1]][3] := 0, !!Z_name[5][[1]][4] := Z5,
  )

n_dim <- 2 + dim_X + dim_Y + dim_Z
end <- dim(df1)[2]
start <- end - 4* n_dim + 1

mlogit.df <- mlogit.data(df1,
                         choice = 'choice',
                         shape = 'wide',
                         varying = start:end,
                         sep="")


model_0 <- mlogit(choice ~ mu + rho  | 0, 
                  data = mlogit.df)
summary(model_0)

model_summary <- summary(model_0)

left_CI <- model_summary$CoefTable[,1]-1.96*model_summary$CoefTable[,2]
right_CI <- model_summary$CoefTable[,1]+1.96*model_summary$CoefTable[,2]
left_CI
right_CI

model_1 <- mlogit(choice ~ mu + rho + age + years_with_firm + status +office+practice+gender+law_school | 0, 
                  data = mlogit.df)
summary(model_1)

model_summary <- summary(model_1)

left_CI <- model_summary$CoefTable[,1]-1.96*model_summary$CoefTable[,2]
right_CI <- model_summary$CoefTable[,1]+1.96*model_summary$CoefTable[,2]
print(left_CI)
print(right_CI)






###### Analysis of trade network data

Log_of_Gravity <- read_excel("trade_network/Log of Gravity.xls")
countrycodes <- read_excel("trade_network/countrycodes.xls")
names(countrycodes)[1] <- "code"

# Significant trade partnerships -------------------------------------------------

# Find signigicant trade partnerships based on trade volume
# Calculate GDP per country
gdps <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lypim)%>%
               group_by(s1_im)%>%
               summarise(log_gdp = mean(lypim))%>%
               mutate(importer_GDP = exp(log_gdp))%>%
               select(s1_im, importer_GDP),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))

Total_import <- Log_of_Gravity%>%
  group_by(s1_im)%>%
  summarise(Total_import = sum(trade))
Total_export <- Log_of_Gravity%>%
  group_by(s2_ex)%>%
  summarise(Total_export = sum(trade))
Total_import_export <- Total_import%>%
  right_join(Total_export, by = c("s1_im" = "s2_ex"))%>%
  right_join(countrycodes, by = c("s1_im" = "code"))%>%
  select(c(1,4,2,3))
names(Total_import_export)[1] <- "code"

Total_import_export$total_volume <- Total_import_export$Total_import + Total_import_export$Total_export
Total_import_export$code <- as.character(Total_import_export$code)

Network_data <- Log_of_Gravity%>%
  select(s1_im, s2_ex, trade)%>%
  right_join(countrycodes[,1:2], by = c("s1_im" = "code"))%>%
  rename(importer_country = country)%>%
  right_join(countrycodes[,1:2], by = c("s2_ex" = "code"))%>%
  rename(exporter_country = country)%>%
  right_join(Total_import, by = "s1_im")%>%
  right_join(Total_export, by = "s2_ex")%>%
  select(s1_im, s2_ex, importer_country, exporter_country,
         trade, Total_import, Total_export)


country_pairs <- Network_data %>%
  rowwise() %>%
  mutate(pair = paste0(sort(c(s1_im, s2_ex)), collapse = "-")) %>%
  ungroup()


total_trade_values <- country_pairs %>%
  group_by(pair) %>%
  summarise(total_trade_value = sum(trade)) %>%
  ungroup() %>%
  separate(pair, into = c("country1", "country2"), sep = "-") %>%
  left_join(Total_import_export, by = c("country1" = "code")) %>%
  rename(country1_total_volume = total_volume)%>%
  left_join(Total_import_export, by = c("country2" = "code")) %>%
  rename(country2_total_volume = total_volume)%>%
  select(country1, country2,total_trade_value,country1_total_volume,country2_total_volume)%>%
  mutate(country1to2 = ifelse(total_trade_value >= 0.01*country1_total_volume, 1, 0),
         country2to1 = ifelse(total_trade_value >= 0.01*country2_total_volume, 1, 0))

link <- total_trade_values %>%
  filter((country1to2 == 1) | (country2to1 == 1))%>%
  select(country1, country2,country1to2,country2to1)%>%
  arrange(country1)

edges <- c()

# Add edges from A to B where A_to_B is 1
edges <- c(edges, as.vector(t(link[link$country1to2 == 1, c("country1", "country2")])))

# Add edges from B to A where B_to_A is 1
edges <- c(edges, as.vector(t(link[link$country2to1 == 1, c("country2", "country1")])))

# Create the graph using the edge list
g <- graph(edges, directed = TRUE)

# Generate the adjacency matrix
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

# Print the adjacency matrix
print(adj_matrix)

print(sprintf("Total number of edges: %s", length(E(g)))) #2141
print(sprintf("Edge density: %s", edge_density(g))) #0.11661220043573

degrees <- rowSums(adj_matrix)


sum(adj_matrix * t(adj_matrix))/2  # 260 

gdps <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lypim)%>%
               group_by(s1_im)%>%
               summarise(log_gdp = mean(lypim))%>%
               mutate(importer_GDP = exp(log_gdp))%>%
               select(s1_im, importer_GDP),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))

remoteness <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lremot_im)%>%
               group_by(s1_im)%>%
               summarise(log_remoteness = mean(lremot_im))%>%
               mutate(importer_remoteness = exp(log_remoteness))%>%
               select(s1_im, importer_remoteness),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))


log_remoteness <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lremot_im)%>%
               group_by(s1_im)%>%
               summarise(log_remoteness = mean(lremot_im))%>%
               select(s1_im, log_remoteness),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))

log_gdps <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, lypim)%>%
               group_by(s1_im)%>%
               summarise(log_gdp = mean(lypim))%>%
               select(s1_im, log_gdp),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))


landlocked <- countrycodes%>%
  right_join(Log_of_Gravity%>%
               select(s1_im, landl_im)%>%
               group_by(s1_im)%>%
               summarise(landlocked = mean(landl_im))%>%
               select(s1_im, landlocked),
             by = c("code" = "s1_im"))%>%
  select(-(3:4))

n <- nrow(countrycodes)
covariates <- Log_of_Gravity%>%select(s1_im, s2_ex, ldist, border, comlang, colony, comfrt)%>%
  filter(FALSE)
for (i in 1:(n-1)) {
  if(i %% 50 == 0){print(i)}
  country_i <- countrycodes$code[i]
  for (j in (i+1):n) {
    country_j <- countrycodes$code[j]
    covariates <- rbind(covariates,
                        Log_of_Gravity%>%select(s1_im, s2_ex, ldist, border, comlang, colony, comfrt)%>%
                          filter(s1_im == country_i, s2_ex == country_j))
  }
}
rm(country_i, country_j, i, j)

node_order <- as.character(landlocked$code)

adj_matrix_reordered <- adj_matrix[node_order, node_order]

node <- dim(adj_matrix_reordered)[1]

dim_X <- 1
dim_Y <- 1
dim_Z <- 5

n_col <- 2 + 2*(dim_X+dim_Y)+ dim_Z
n_row <- node*(node-1)/2

data <- matrix(rep(0,n_row*n_col),nrow = n_row)

for (i in c(1:(node-1))) {
  for (j in c((i+1):node)){
    A_ij <- adj_matrix_reordered[i,j]
    A_ji <- adj_matrix_reordered[j,i]
    X_i <- as.vector(unlist(landlocked[i,3]))
    X_j <- as.vector(unlist(landlocked[j,3]))
    Y_i <- as.vector(unlist(log_gdps[i,3]))
    Y_j <- as.vector(unlist(log_gdps[j,3]))
    edge_index <- (j-i)+(i-1)*(2*node-i)/2
    Z_ij <- as.vector(unlist(covariates[edge_index,3:7]))
    data[edge_index,] <- c(A_ij,A_ji,X_i,X_j,Y_i,Y_j,Z_ij)
  }
}


df1 <- as.data.frame(data) %>%
  rename(
    A_ij = V1, A_ji = V2,
    X_i = V3, X_j = V4,
    Y_i = V5, Y_j = V6,
    Z1 = V7, Z2 = V8, Z3 = V9, Z4 = V10, Z5 = V11
  )


X_name <- rep(0,4)
Y_name <- rep(0,4)
Z_name <- list()
for (i in 1:4){
  X_name[i] <- paste0('landlocked',i)
  Y_name[i] <- paste0('gdp',i)
}
for (j in 1:5){
  Z <- rep(0,4)
  for (i in 1:4){
    Z[i] <- paste0(names(covariates)[2+j],i)
  }
  Z_name <- append(Z_name, list(Z))
}


df1 <- df1 %>%
  mutate(
    choice = case_when(
      A_ij == 0 & A_ji == 0 ~ 1,
      A_ij == 1 & A_ji == 0 ~ 2,
      A_ij == 0 & A_ji == 1 ~ 3,
      TRUE ~ 4
    ),
    mu1 = 0, mu2 = 1, mu3 = 1, mu4 = 2,
    rho1 = 0, rho2 = 0, rho3 = 0, rho4 = 1,
    !!X_name[1] := 0,!!X_name[2] := X_i, !!X_name[3] := X_j, !!X_name[4] := X_i + X_j,
    !!Y_name[1] := 0, !!Y_name[2] := Y_j, !!Y_name[3] := Y_i, !!Y_name[4] := Y_i + Y_j,
    !!Z_name[1][[1]][1] := 0, !!Z_name[1][[1]][2] := 0, !!Z_name[1][[1]][3] := 0, !!Z_name[1][[1]][4] := Z1,
    !!Z_name[2][[1]][1] := 0, !!Z_name[2][[1]][2] := 0, !!Z_name[2][[1]][3] := 0, !!Z_name[2][[1]][4] := Z2,
    !!Z_name[3][[1]][1] := 0, !!Z_name[3][[1]][2] := 0, !!Z_name[3][[1]][3] := 0, !!Z_name[3][[1]][4] := Z3,
    !!Z_name[4][[1]][1] := 0, !!Z_name[4][[1]][2] := 0, !!Z_name[4][[1]][3] := 0, !!Z_name[4][[1]][4] := Z4,
    !!Z_name[5][[1]][1] := 0, !!Z_name[5][[1]][2] := 0, !!Z_name[5][[1]][3] := 0, !!Z_name[5][[1]][4] := Z5,
  )

n_dim <- 2 + dim_X + dim_Y + dim_Z
end <- dim(df1)[2]
start <- end - 4* n_dim + 1

mlogit.df <- mlogit.data(df1,
                         choice = 'choice',
                         shape = 'wide',
                         varying = start:end,
                         sep="")

model_0 <- mlogit(choice ~ mu + rho  | 0, 
                  data = mlogit.df)
summary(model_0)

model <- mlogit(choice ~ mu + rho +gdp+landlocked+ ldist+border +comlang +colony+comfrt| 0, 
                data = mlogit.df)

summary(model)

model_summary <- summary(model)

left_CI <- model_summary$CoefTable[,1]-1.96*model_summary$CoefTable[,2]
right_CI <- model_summary$CoefTable[,1]+1.96*model_summary$CoefTable[,2]

print(left_CI)
print(right_CI)






