
#################
##### SETUP #####
#################

# cleaning the environment
rm(list = ls())

# Uploading/Installing the packages
list.of.packages <- c("readr", "tidyverse", "regclass", "huxtable", "jtools",
                      "stats","olsrr", "glmnet")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

invisible(lapply(list.of.packages, library, character.only = TRUE))

# Getting the dataset (from my public repository)
data="https://raw.githubusercontent.com/carlosleong/climate_change_test/main/climate_change.csv"
data<-read_csv(url(data))

##################################
######## Data Managment ##########
##################################

#Changing var names patterns (bcs problematic)
for (i in 1:length(names(data))){
aux <- names(data[,i])
names(data)[i] <- gsub("-", "_", aux)
}

#Defining set of covariates and independent variables
Covariates <- setdiff(names(data), c("Temp", "Year", "Month"))

#splitting into training and test datasets

data_train <- data[data$Year<=2006,] 
data_test <-  data[data$Year>2006,] 

##################################
########  Analysis ###############
##################################

#Linea Model
f1 <- as.formula(paste("Temp", paste(Covariates, 
                  collapse = ' + '), sep = " ~ "))

linear_1 <- lm(f1, data_train)
summary(linear_1)


#Correlations among Covariates and N2O & CFC_11 
aux <- c("N2O", "CFC_11") 
corr_matrix <- cor(data_train) %>%
               subset(,c(aux,setdiff(Covariates,aux))) 

corr_matrix <- corr_matrix[aux,]
corr_matrix

# Common trend
data_train$year_f <- as.factor(data_train$Year)
data_train$month_f <- as.factor(data_train$Month)
aux <- c("Temp", "N2O", "CFC_11")

for (i in 1:length(aux)){
 assign(paste("lm_aux",i,sep = ""),
        lm(get(aux[i]) ~ year_f , data = data_train))  
}
export_summs(lm_aux1, lm_aux2, lm_aux3, model.names = c("Temp","N2O","CFC_11"))

###############################################################
################### Simplifying models ########################
###############################################################
# Matrix of best specs (training and test data)
WINNERS <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 5))
colnames(WINNERS) <- c("Specification", "AIC-ante", "R2-ante", "R2adj-ante",
                       "R2-post")
rownames(WINNERS) <- c("Best Subset", "Forward Stepwi", "Ridge", "Lasso")

# Functions that will be used to assess the predictive power in the test data, and in the case of
# the glmnet training data too (bcs the program don't perform the statistic)

#R2
R2 <- function(true, predicted){
  SSE <- sum((predicted-true)^2)
  SST <- sum((true-mean(true))^2)
  R2 <- 1 - (SSE/SST)
  return(R2)
}

#R2adj
R2adj <- function(R2, N, k){
  R2adj <- 1 - ((1-R2)*(N-1)/(N-k-1))
  return(R2adj)
} 

#R2-predicted
predict_R2 <- function(true_train, true_test, predicted){
  SSE <- sum((predicted-true_test)^2)
  SST <- sum((true_test-mean(true_train))^2)
  predict_R2 <- 1 - (SSE/SST)
  return(predict_R2)
}


#########
#Method 1: Best subset selection 
#########

# This section implements the algorithm proposed by:
# James, G.; Witten, D.; Hastie, T.; and Tibshirani, R. (2013) in their textbook
# An Introduction to Statistical Learning, Edited by Springers

# Algorithm 6.1 (pp. 205)
model_winnerset <- as.data.frame(matrix(data = 0, nrow = length(Covariates)+1, ncol = 4))
colnames(model_winnerset) <- c("Specification", "AIC","R2", "R2adj")

#Step 1
model_winnerset[1,1] <- "Temp~1"
aux <- lm(Temp~1, data_train)    
model_winnerset[1,2] <- AIC(aux)
aux <- summary(aux)
model_winnerset[1,3] <- aux[["r.squared"]]
model_winnerset[1,4] <- aux[["adj.r.squared"]]

#Step 2
for (i in 1:length(Covariates)){
  aux <- combn(1:length(Covariates),i, simplify = TRUE)
  model_contestant <- as.data.frame(matrix(data = 0, nrow = ncol(aux), ncol = 4))
for (j in 1:ncol(aux)){
    model <- aux[,j]
    model <- Covariates[model]
    model <- paste("Temp", paste(model,collapse = ' + '), sep = " ~ ")
    model_contestant[j,1] <- model
    model <- lm(as.formula(model), data_train)
    model_contestant[j,2] <- AIC(model)
    model <- summary(model)
    model_contestant[j,3] <- model[["r.squared"]]
    model_contestant[j,4] <- model[["adj.r.squared"]]
}
    max_r2 <- max(model_contestant[,3])
    max_r2 <- grep(max_r2, model_contestant[,3])
    model_winnerset[i+1,] <-  model_contestant[max_r2,]
}

#Step 3
    min_AIC <- min(model_winnerset[,2])
    min_AIC <- grep(min_AIC, model_winnerset[,2])
    WINNERS[1,1] <- model_winnerset[min_AIC,1]
    WINNERS[1,2] <- model_winnerset[min_AIC,2]
    WINNERS[1,3] <- model_winnerset[min_AIC,3]
    WINNERS[1,4] <- model_winnerset[min_AIC,4]
    
#Prediction (test data)
best_subs_model <- lm(as.formula(WINNERS[1,1]),data = data_train)
best_subs_model_predict <- predict(best_subs_model, newdata=data_test)

WINNERS[1,5] <- predict_R2(data_train$Temp,data_test$Temp, best_subs_model_predict)    

  
#########
#Method 2: Forward stepwise selection 
#########
fwd_stepmod  <- lm(f1, data = data_train)
fwd_stepmod  <- ols_step_forward_p(fwd_stepmod)
fwd_stepmod  <- fwd_stepmod[["predictors"]]
fwd_stepmod  <- paste("Temp", paste(fwd_stepmod,collapse = ' + '), sep = " ~ ")
WINNERS[2,1] <- fwd_stepmod
fwd_stepmod  <- lm(as.formula(fwd_stepmod), data = data_train)
WINNERS[2,2] <- AIC(fwd_stepmod)
fwd_stepmod  <- summary(fwd_stepmod)    
WINNERS[2,3] <- fwd_stepmod[["r.squared"]] 
WINNERS[2,4] <- fwd_stepmod[["adj.r.squared"]]   
   

#Prediction (test data)
fwd_stepmod <- lm(as.formula(WINNERS[2,1]),data = data_train)
fwd_stepmod_predict <- predict(fwd_stepmod, newdata=data_test)

WINNERS[2,5] <- predict_R2(data_train$Temp,data_test$Temp, fwd_stepmod_predict)    

    
#########
#Method 3: Ridge Regression 
#########

# Step 1: find optimal lambda
ridge_model <- glmnet(as.matrix(data_train[Covariates]), data_train$Temp,
                      family="gaussian", alpha = 0)

opt_lambda_ridge <- min(ridge_model[["lambda"]])

# Step 2: Estimate the training model
predict_ridge_train <- predict(ridge_model, s = opt_lambda_ridge, 
                             newx = as.matrix(data_train[Covariates]))

# Storing results
WINNERS[3,1] <- "All covariates included"
WINNERS[3,2] <- "NA"
WINNERS[3,3] <- R2(data_train$Temp, predict_ridge_train)
WINNERS[3,4] <- R2adj(WINNERS[3,3], nrow(data_train), length(Covariates)+1)

#Prediction (test data)
predict_ridge_test <- predict(ridge_model, s = opt_lambda_ridge, 
                              newx = as.matrix(data_test[Covariates]))

WINNERS[3,5] <- predict_R2(data_train$Temp,data_test$Temp, predict_ridge_test)    


#########
#Method 4: Lasso Regression 
#########

# Step 1: find optimal lambda
lasso_model <- glmnet(as.matrix(data_train[Covariates]), data_train$Temp,
                      family="gaussian", alpha = 1)
opt_lambda_lasso <- min(lasso_model[["lambda"]])

# Step 2: Estimate the training model
predict_lasso_train <- predict(lasso_model, s = opt_lambda_lasso, 
                              newx = as.matrix(data_train[Covariates]))

# Storing results
WINNERS[4,1] <- "All covariates included"
WINNERS[4,2] <- "NA"
WINNERS[4,3] <- R2(data_train$Temp, predict_lasso_train)
WINNERS[4,4] <- R2adj(WINNERS[3,3], nrow(data_train), length(Covariates)+1)

#Prediction (test data)
predict_lasso_test <- predict(lasso_model, s = opt_lambda_lasso, 
                              newx = as.matrix(data_test[Covariates]))

WINNERS[4,5] <- predict_R2(data_train$Temp,data_test$Temp, predict_lasso_test)    

WINNERS[1:4,]



