# dealing with imbalanced data
library(ROCR)
library(grid)
library(broom)
library(caret)
library(tidyr)
library(dplyr)
library(scales)
library(ggplot2)
library(ggthemr) 
library(ggthemes)
library(gridExtra)
library(data.table)
# based on previous 
setwd('C:/Users/46541/Desktop/CEE Research/old_data') 
data <- read.csv("ROC_data.csv")  # very large!
head(data)
res_col <- c(4, 5)
previous_result <- data[,res_col]

mydata <- data[,-res_col]
head(mydata)
dim(mydata)
summary(mydata)
# DATA CLASSES
# LS Types: 8 = No landslides, 1 = Fall/Topple, 2 = Torrent, 
#       3 = Avalanche, 4 = Slump/Creep, 5 = Sackung
# Debris avalanche areas = 1 runout, 2 source, 3 other
# Lithology: 1 = Unconsolidated Sediment, 2 = Ultrabasic rock,
#  3 = Weak Metamorphic Foliated, 4 = Sedimentary Rock,
#  5 = Hard Metamorphic, 6 = Intrusive Igneous, 
#  7 = Volcanic/Extrusive Igneous
# LULC types: 71 = Herbaceous, 52 = Shrubland, 41 = Forest, 31 = Barren,
#  21 = Developed

# cleaning data
sum(is.na(mydata))

# transform LS_type into binary class : 
# 0: no landslide happen
# 1: landslide happens
mydata$LS_type[mydata$LS_type!=8] <- 1
mydata$LS_type[mydata$LS_type==8] <- 0
mydata$LS_type <- as.factor(mydata$LS_type)

table(mydata$LS_type)
prop.table(table(mydata$LS_type))

#--------------------------------------------------------------------------------------------------
# split the whole data into 2 pieces: 80% for training, and 20% for testing 
# we will use training data to trai the model
# and use the same testing data to find which model performs the best
set.seed(123)
test_index <- createDataPartition(mydata$LS_type, p= 0.2, list =FALSE)
# training data
data_train <- mydata[-test_index, ]
table(data_train$LS_type)
prop.table(table(data_train$LS_type))
# testing data
data_test <- mydata[test_index, ]
table(data_test$LS_type)
prop.table(table(data_test$LS_type))
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
# Baselines of our models
#      the baseline is always not satisfied, and our goal is to improve it!
#--------------------------------------------------------------------------------------------------

# Model1: logistic regression
train_logit <- data_train[,-c(1,3)]
test_logit <- data_test[, -c(1,3)]

train_logit_baseline <- train_logit 
test_logit_baseline <- test_logit
str(train_logit_baseline)
# we need to transform litho and lulc into factors
train_logit_baseline$litho <- as.factor(train_logit_baseline$litho)
train_logit_baseline$lulc <- as.factor(train_logit_baseline$lulc)
test_logit_baseline$litho <- as.factor(test_logit_baseline$litho)
test_logit_baseline$lulc <- as.factor(test_logit_baseline$lulc)

glm_baseline <-  glm(LS_type~., data = train_logit_baseline, family= binomial)
summary_glm_baseline = summary(glm_baseline )
list( summary_glm_baseline$coefficient, 
      round( 1 - ( summary_glm_baseline$deviance / summary_glm_baseline$null.deviance ),3) )
pred_baseline <- predict(glm_baseline, newdata = test_logit_baseline, type = "response")

# baseline results
#-----------------------------------------------------------
# confusion_info is a function returns some useful metrices 
confusion_info <- function(prediction, labels, cutoff=0.5){
  index <- which(prediction == 0)
  if (length(index) == 0){
    if (sum(prediction>cutoff) == 0){
      return(paste("There are no predictions > cutoff when cutoff =", cutoff))
    }
    pred.class <- as.integer(prediction  > cutoff)
    cft <- table(pred.class, labels)
    tp <- cft[2, 2]
    tn <- cft[1, 1]
    fp <- cft[2, 1]
    fn <- cft[1, 2]
    accuracy <- (tp + tn)/(tp + tn + fp + fn)
    sensitivity <- tp/(tp + fn)  # recall
    recall_val <- sensitivity
    precision_val <- tp/(tp + fp) # precision
    specificity <- tn/(tn + fp)
    f1_score <- 2*precision_val*recall_val/(precision_val+recall_val)
    
    return(list(confusion_matrix_table = cft,
                Accuracy = accuracy,
                precision = precision_val,
                Sensitivity = sensitivity,
                Specificity = specificity,
                f1_score = f1_score
    ))
    
  }
  else{
    cft <- table(prediction, labels)
    tp <- cft[2, 2]
    tn <- cft[1, 1]
    fp <- cft[2, 1]
    fn <- cft[1, 2]
    accuracy <- (tp + tn)/(tp + tn + fp + fn)
    sensitivity <- tp/(tp + fn)  # recall
    recall_val <- sensitivity
    precision_val <- tp/(tp + fp) # precision
    specificity <- tn/(tn + fp)
    f1_score <- 2*precision_val*recall_val/(precision_val+recall_val)
    
    return(list(confusion_matrix_table = cft,
                Accuracy = accuracy,
                Sensitivity = sensitivity,
                precision = precision_val,
                Specificity = specificity,
                f1_score = f1_score
    ))
  }
}
#----------------------------------------------------------------
confusion_info(pred_baseline, test_logit_baseline$LS_type, cutoff = 0.01)

#----------------------------------------------------------------------------------------
# plots and results of baseline
cutoff <- seq(0.01,max(pred_baseline),0.01)
accuracy_baseline <- cutoff
f1_score_baseline <- cutoff
confusion_info_accuracy <- function(cutoff){
  c1 <- confusion_info(pred_baseline,test_logit_baseline$LS_type,cutoff)
  return(c1$Accuracy)
}
accuracy_baseline <- sapply(cutoff, confusion_info_accuracy)
acb <- data.frame(x = cutoff, y = accuracy_baseline)
ggplot(acb, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='Accuracy',title='Baseline accuracy (Logistic regression)')
  
confusion_info_f1 <- function(cutoff){
  c1 <- confusion_info(pred_baseline,test_logit_baseline$LS_type,cutoff)
  return(c1$f1_score)
}
f1_score_baseline <- sapply(cutoff, confusion_info_f1)
f1b <- data.frame(x = cutoff[1:30], y = f1_score_baseline[1:30])
ggplot(f1b, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='Baseline f1 score (Logistic regression)')

paste("The optimal accuracy is", max(accuracy_baseline), 
      "and the corresponding cutoff value is", cutoff[accuracy_baseline == max(accuracy_baseline)][1] )
confusion_info(pred_baseline, test_logit_baseline$LS_type, 
               cutoff = cutoff[accuracy_baseline == max(accuracy_baseline)][1])

paste("The optimal f1 score is",max(f1_score_baseline[1:30]), 
      "and the corresponding cutoff value is",cutoff[8] )
confusion_info(pred_baseline, test_logit_baseline$LS_type, cutoff = cutoff[8])

# ROC AUC P-R curve
library(ROSE)
roc.curve(test_logit_baseline$LS_type, pred_baseline)


#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
# Model 1*: logistic regression with Lasso, Ridge
# in order to avoid overfitting
train_logit_lasso <- train_logit
test_logit_lasso <- test_logit

train_logit_lasso$litho <- as.factor(train_logit_lasso$litho)
train_logit_lasso$lulc <- as.factor(train_logit_lasso$lulc)
xfactors<-model.matrix(train_logit_lasso$LS_type~train_logit_lasso$litho+train_logit_lasso$lulc)[,-1]
x<-as.matrix(data.frame(train_logit_lasso[,c(2,3,4,5,6)],xfactors))
glmmod<-glmnet(x,y=train_logit_lasso$LS_type, alpha=1, family='binomial', lambda = seq(0,1,0.00001))
glmmod
plot(glmmod,xvar="lambda")
coef(glmmod)[,10]
#-----------------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------------
# Model 2.1: Random Forest
library(randomForest)
train_RF <- train_logit
test_RF <- test_logit
str(train_RF)
# ensure that the response is factor data, predictors are numeric
# we need to under sample first
under_sample_merge <- function(data, multiply_1=1){
  # data: raw data
  # multiply_1 : the ratio of multiple class to minority class we set
  d_1 <- data[data$LS_type == 1,]
  num_1 <- dim(d_1)[1]
  d_remaining <- data[data$LS_type != 1,]
  num_0 <- dim(d_remaining)[1]
  pick_index <- sample(num_0, multiply_1*num_1)
  d_pick <- d_remaining[pick_index,]
  # smaller new dataset
  d_merge <- rbind(d_1, d_pick)
  return(d_merge)
}
set.seed(1)
train_under_RF <- under_sample_merge(train_RF, 1)
dim(train_under_RF)
#RF_baseline <- randomForest(LS_type~., data=train_under_RF_reduce, importance=TRUE)
RF_baseline <- randomForest(y=train_under_RF$LS_type, x=train_under_RF[,-1], data=train_under_RF)

pre_RF_baseline <- predict(RF_baseline, test_RF)

# result
confusion_info(pre_RF_baseline, test_RF$LS_type)

# optimization
print(RF_baseline)
plot(RF_baseline)

# to find the optimal nodes
rate <- seq(1,ncol(train_under_RF_reduce)/2)
for(i in 1:ncol(train_under_RF_reduce)/2){
  set.seed(123)
  model.forest <-randomForest(LS_type ~ ., data = train_under_RF,mtry=i,importance=TRUE,ntree=400)
  rate[i] = mean(model.forest$err.rate[300,])
}

set.seed(1)
model.forest1 <- randomForest(y=train_under_RF$LS_type, x=train_under_RF[,-1], data=train_under_RF,mtry = 2,ntree=500)
pre_RF_baseline <- predict(model.forest1, test_RF)
confusion_info(pre_RF_baseline, test_RF$LS_type)

# 1:0.3264 2:0.3309 3:0.3280 4:0.33617 5:0.33211 6: 0.3337
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
# Model 2.2: Adaboost

#-----------------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------------------
# resample methods
#-----------------------------------------------------------------------------------------------------
# transform data type in test_logit
str(test_logit)
test_logit$litho <- as.factor(test_logit$litho)
test_logit$lulc <- as.factor(test_logit$lulc)
str(test_logit)
#-----------------------------------------------------------------------------------------------------
# resample methods1: random over sample
library(ROSE)
data_balanced_over <- ovun.sample(LS_type ~ ., data = train_logit, method = "over",N = 3000000)$data
table(data_balanced_over$LS_type)
prop.table(table(data_balanced_over$LS_type))
# as.factor litho and lulc
data_balanced_over$litho <- as.factor(data_balanced_over$litho)
data_balanced_over$lulc <- as.factor(data_balanced_over$lulc)
model_glm_over <- glm( LS_type ~. , data_balanced_over, family = binomial)
pred_glm_over <- predict(model_glm_over, newdata = test_logit, type = 'response')

#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_over, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.2,0.8,0.01)
f1_over <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_over)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='Random oversampling')
cutoff_opt <- cutoff[which(f1_over == max(f1_over))]
paste("The optimal f1 score is",max(f1_over),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_over, test_logit$LS_type,cutoff=cutoff_opt)

library(pROC)
roc_over <- roc(test_logit$LS_type, pred_glm_over)
auc(roc_over)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Randomly under sampling   
data_balanced_under <- ovun.sample(LS_type ~ ., data = train_logit, method = "under", N = 140000, seed = 1)$data
table(data_balanced_under$LS_type)
prop.table(table(data_balanced_under$LS_type))
# as.factor
data_balanced_under$litho <- as.factor(data_balanced_under$litho)
data_balanced_under$lulc <- as.factor(data_balanced_under$lulc)
str(data_balanced_under)
model_under <- glm( LS_type ~. , data_balanced_under, family = binomial)

pred_glm_under <- predict(model_under, newdata = test_logit, type = 'response')
#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_under, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.2,0.8,0.01)
f1_under <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_under)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='Random oversampling')
cutoff_opt <- cutoff[which(f1_under == max(f1_under))]
paste("The optimal f1 score is",max(f1_over),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_under, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_glm_under)
auc(roc_under)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Under-sampling + Over-sampling
data_balanced_both <- ovun.sample(LS_type ~ ., data =train_logit, method = "both", p=0.45, N=250000, seed = 1)$data
table(data_balanced_both$LS_type)
prop.table(table(data_balanced_both$LS_type))
# as.factor
data_balanced_both$litho <- as.factor(data_balanced_both$litho)
data_balanced_both$lulc <- as.factor(data_balanced_both$lulc)
model_both <- glm( LS_type ~. , data_balanced_both, family = binomial)
pred_glm_both <- predict(model_both, newdata = test_logit, type = 'response')
#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_both, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.2,0.8,0.01)
f1_both <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_both)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='Both')
cutoff_opt <- cutoff[which(f1_both == max(f1_both))]
paste("The optimal f1 score is",max(f1_both),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_both, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_glm_both)
auc(roc_under)

#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# ROSE
data_balanced_rose <- train_logit
str(data_balanced_rose)
data_balanced_rose$litho <- as.factor(data_balanced_rose$litho )
data_balanced_rose$lulc <- as.factor(data_balanced_rose$lulc)

str(data_balanced_rose)
data_balanced_rose <- ROSE(LS_type ~ ., data = data_balanced_rose , seed = 1)$data
str(data_balanced_rose)

model_glm_rose <- glm( LS_type ~. , data_balanced_rose, family = binomial)
pred_glm_rose <- predict(model_glm_rose, newdata = test_logit, type = 'response')
#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_rose, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.1,0.9,0.01)
f1_both <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_both)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='ROSE')
cutoff_opt <- cutoff[which(f1_both == max(f1_both))]
paste("The optimal f1 score is",max(f1_both),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_rose, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_glm_rose)
auc(roc_under)
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
# Tomek link
library(unbalanced)
output <-train_logit$LS_type
input <- train_logit[,-1]
str(input)
data_tlink <- ubTomek(X=input, Y=output)
data_tlink <- cbind(data_tlink$X,data_tlink$Y)
colnames(data_tlink)[8] <- 'LS_type'
table(data_tlink$LS_type)
prop.table(table(data_tlink$LS_type))
# as.factor
data_tlink$lulc <- as.factor(data_tlink$lulc)
data_tlink$litho <- as.factor(data_tlink$litho)

model_glm_tlink <- glm( LS_type ~ . , data = data_tlink, family = binomial )
pred_glm_tlink <- predict(model_glm_tlink, newdata = test_logit, type = 'response')

#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_tlink, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.001,0.3,0.01)
f1 <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='Tomek link')
cutoff_opt <- cutoff[which(f1 == max(f1))]
paste("The optimal f1 score is",max(f1),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_tlink, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_glm_tlink)
auc(roc_under)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# SMOTE
library(DMwR)
set.seed(1)
train_under <- under_sample_merge(train_logit, 10)
table(train_under$LS_type)
prop.table(table(train_under$LS_type))

#train_random <- train_under[sample(dim(train_under)[1],10000),]
#table(train_random$LS_type)
#prop.table(table(train_random$LS_type))
train_under$litho <- as.factor(train_under$litho)
train_under$lulc <- as.factor(train_under$lulc)

data_SMOTE <- SMOTE(LS_type ~., train_under, perc.over = 450,  perc.under = 250)
table(data_SMOTE$LS_type)
prop.table((table(data_SMOTE$LS_type)))

dim(data_SMOTE)
# since we might apply SMOTE + Tomek link, we will not as.factor in data_SMOTE

model_glm_smote <- glm( LS_type ~ . , data = data_SMOTE, family = binomial )
pred_glm_smote <- predict(model_glm_smote, newdata = test_logit, type = 'response')

#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_smote, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.1,0.8,0.01)
f1 <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='SMOTE')
cutoff_opt <- cutoff[which(f1 == max(f1))]
paste("The optimal f1 score is",max(f1),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_smote, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_glm_smote)
auc(roc_under)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# SMOTE + T-link
str(data_SMOTE)

output <- data_SMOTE$LS_type
input <- data_SMOTE[,-1]
input$litho <- as.numeric(input$litho)
input$lulc <- as.numeric(input$lulc)
str(input)
data_smote_tlink <- ubTomek(X=input, Y=output)
data_smote_tlink <- cbind(data_smote_tlink$X,data_smote_tlink$Y)
colnames(data_smote_tlink)[8] <- 'LS_type'
table(data_smote_tlink$LS_type)
prop.table(table(data_smote_tlink$LS_type))
data_smote_tlink$litho <- as.factor(data_smote_tlink$litho)
data_smote_tlink$lulc <- as.factor(data_smote_tlink$lulc)

model_glm_smotetlink <- glm( LS_type ~ . , data = data_smote_tlink, family = binomial )

str(test_logit)
test_smote_family <- test_logit
test_smote_family$litho <- as.numeric(test_smote_family$litho)
test_smote_family$lulc <- as.numeric(test_smote_family$lulc)
test_smote_family$litho <- as.factor(test_smote_family$litho)
test_smote_family$lulc <- as.factor(test_smote_family$lulc)
str(test_smote_family)

pred_glm_smotetlink <- predict(model_glm_smotetlink, newdata = test_smote_family, type = 'response')

#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_glm_smotetlink, test_smote_family$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.2,0.8,0.01)
f1_both <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_both)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='SMOTE+T-link')
cutoff_opt <- cutoff[which(f1_both == max(f1_both))]
paste("The optimal f1 score is",max(f1_both),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_glm_smotetlink, test_smote_family$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_smote_family$LS_type, pred_glm_smotetlink)
auc(roc_under)

#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Rndom Undersampling + T-link
data_under_tlink <- ovun.sample(LS_type ~ ., data = train_logit, method = "under", N = 250000, seed = 123)$data
table(data_under_tlink$LS_type)

output <-data_under_tlink$LS_type
input <- data_under_tlink[,-1]
str(input)
data_under_tlink <- ubTomek(X=input, Y=output)
data_under_tlink <- cbind(data_under_tlink$X,data_under_tlink$Y)
colnames(data_under_tlink)[8] <- 'LS_type'

data_under_tlink$litho <- as.factor(data_under_tlink$litho)
data_under_tlink$lulc <- as.factor(data_under_tlink$lulc)
str(data_under_tlink)
model_under_tlink <- glm(LS_type~., data = data_under_tlink, family = binomial)
table(data_under_tlink$LS_type)
prop.table(table(data_under_tlink$LS_type))
pred_under_tlink <- predict(model_under_tlink, newdata = test_logit, type = 'response')
confusion_info(pred_under_tlink, test_logit$LS_type)

#plots
c_f1 <- function(cutoff){
  f1 <- confusion_info(pred_under_tlink, test_logit$LS_type,cutoff)$f1_score
}
cutoff <- seq(0.1,0.6,0.01)
f1_both <- sapply(cutoff,c_f1)
d <- data.frame(x=cutoff,y=f1_both)
ggplot(d, aes(x=x,y=y))+
  geom_point()+
  geom_line()+
  labs(x='cutoff values', y='f1 score',title='T-link+Undersampling')
cutoff_opt <- cutoff[which(f1_both == max(f1_both))]
paste("The optimal f1 score is",max(f1_both),"and the corresponding cutoff is", cutoff_opt)
confusion_info(pred_under_tlink, test_logit$LS_type,cutoff=cutoff_opt)
roc_under <- roc(test_logit$LS_type, pred_under_tlink)
auc(roc_under)
#-----------------------------------------------------------------------------------------------------







# -------------------------------------------------------------------------
#			Sampling-based method 1:  randomly under-sampling
# step 1: extract from the dataset all the cases labelled as "1" (i.e. landslide happen)
# step 2: randomly select an equal number of cases labelled as "0" from the 
# remaining data
# step 3: merge them into a smaller dataset, but a balanced one
# -------------------------------------------------------------------------
data_pick_merge_partition <- function(data, multiply_1=1 ,ratio_partition=0.25){
  # data: raw data
  # multiply_1 : the nymber we want to enlarge the "0" class
  # ratio_partition: the ratio we want to divide the data into training and validation set
  d_1 <- data[data$LS_type == 1,]
  num_1 <- dim(d_1)[1]
  d_remaining <- data[data$LS_type != 1,]
  num_0 <- dim(d_remaining)[1]
  pick_index <- sample(num_0, multiply_1*num_1)
  d_pick <- d_remaining[pick_index,]
  # smaller new dataset
  d_merge <- rbind(d_1, d_pick)
  # then partition
  id <- createDataPartition(d_merge$LS_type, p =ratio_partition, list = FALSE )
  d_train <- d_merge[-id,]
  d_vali <- d_merge[id,]
  
  return(list(trainig = d_train, 
              validation = d_vali))
}

# note: we here use all of the data
dd <- data_pick_merge_partition(mydata_use, multiply_1=1, ratio_partition=0.25)
d_train <- dd$trainig
prop.table(table(d_train$LS_type) )
d_vali <- dd$validation
prop.table(table(d_vali$LS_type) )
# --------------------------------------------------------
model_glm <- glm( LS_type ~. , d_train, family = binomial )
summary_glm <- summary(model_glm)
summary_glm

# list coefficient and R2
list( summary_glm$coefficient, 
      round( 1 - ( summary_glm$deviance / summary_glm$null.deviance ),3) )

d_train$prediction <- predict( model_glm, newdata = d_train, type = "response" )
d_vali$prediction  <- predict( model_glm, newdata = d_vali , type = "response" )

pred_rand_under <- predict(model_glm, newdata = data_test, type = "response")
confusion_info(pred_rand_under, data_test$LS_type, cutoff=0.5)
roc.curve(data_test$LS_type, pred_rand_under)


confusion_info(d_train$prediction, d_train$LS_type, cutoff=0.5)
confusion_info(d_vali$prediction, d_vali$LS_type, cutoff=0.5)
# without stepAIC or Lasso
library(pROC)
par(pty="s")
roc( d_train$LS_type,d_train$prediction,plot=TRUE ,legacy.axes=TRUE)
roc( d_vali$LS_type,d_vali$prediction,plot=TRUE ,legacy.axes=TRUE)






# -------------------------------------------------------------------------
#			Sampling-based method 2: Tomek link
# -------------------------------------------------------------------------

library(unbalanced)
output2 <-data_train$LS_type
input2 <- data_train[,-1]
input2$litho <- as.numeric(input2$litho )
input2$lulc <- as.numeric(input2$lulc)
dd <- ubTomek(X=input2, Y=output2)
newdd <- cbind(dd$X,dd$Y)
colnames(newdd)[8] <- 'LS_type'
newdd$lulc <- as.factor(newdd$lulc)
newdd$litho <- as.factor(newdd$litho)

model_glm2 <- glm( LS_type ~ . , data = newdd, family = binomial )
summary_glm2 <- summary(model_glm2)
summary_glm2

# list coefficient and R2
list( summary_glm2$coefficient, 
      round( 1 - ( summary_glm2$deviance / summary_glm2$null.deviance ),3) )
# R^2 = 0.058 slightly improved
data_train$prediction <- predict( model_glm2, newdata = data_train, type = "response" )
data_test$prediction  <- predict( model_glm2, newdata = data_test , type = "response" )
library(pROC)
par(pty="s")
roc( data_train$LS_type,data_train$prediction,plot=TRUE ,legacy.axes=TRUE)
roc( data_test$LS_type,data_test$prediction,plot=TRUE ,legacy.axes=TRUE)
# roc_info <- roc( data_train$LS_type,data_train$prediction)




# ------------try a little part of data
data_try <- data[9000:14000,]
data_try$LS_type[data_try$LS_type != 8] <- 1
data_try$LS_type[data_try$LS_type == 8] <- 0
table(data_try$LS_type)
prop.table(table(data_try$LS_type)) 
data_try$X <- c()
data_try$prob_inh <- c()
data_try$prob_dyn <- c()
data_try$DA_area <- c()

library(unbalanced)
outp <-data_try$LS_type
inp <- data_try[,-1]
dddd <- ubTomek(X=inp, Y=outp)
newdddd <- cbind(dddd$X,dddd$Y)
colnames(newdddd)[8] <- 'LS_type'
table(newdddd$LS_type)
prop.table(table(newdddd$LS_type))

library(DMwR)
data_try$LS_type <- as.factor(data_try$LS_type)
balanced_try <- SMOTE(LS_type ~., data_try, perc.over = 1000,  perc.under = 150)
table(balanced_try$LS_type)
# ??????tomek
outp2 <- balanced_try$LS_type
inp2 <- balanced_try[,-1]
balance_reduce <- ubTomek(X=inp2, Y=outp2)
balance_reduce_data <- cbind(balance_reduce$X,balance_reduce$Y )
colnames(balance_reduce_data)[8] <- 'LS_type'
table(balance_reduce_data$LS_type)




#---------------------------

# -------------------------------------------------------------------------
#			Sampling-based method 3: random over-sampling

# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
#			Sampling-based method 4: Smote : Synthetic Minority Oversampling Technique 
# e.g. https://rpubs.com/abhaypadda/smote-for-imbalanced-data
# -------------------------------------------------------------------------
data_train$prediction <- c()
library(DMwR)
balanced.data <- SMOTE(LS_type ~., data_train, perc.over = 100, k = 5, perc.under = 200)


data(iris)
data_eg <- iris[, c(1, 2, 5)]
data_eg$Species <- factor(ifelse(data_eg$Species == "setosa","rare","common")) 
## checking the class distribution of this artificial data set
data_eg$Species[1:40] <- 'common'
table(data_eg$Species)

## now using SMOTE to create a more "balanced problem"
newData <- SMOTE(Species ~ ., data_eg, perc.over = 600,perc.under=200)
table(newData$Species)










