# Final model

##############################################################
# packages
import numpy as np
import pandas as pd
from pandas import Series
from patsy import dmatrices 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler, MaxAbsScaler, MinMaxScaler
from collections import Counter
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score,confusion_matrix
from sklearn.linear_model import LogisticRegression
from imblearn.under_sampling import TomekLinks

##############################################################


##############################################################
### DATA PROCESSSING

#Transforming continuous data into categorical data by CART technique
# functions
def calc_score_median(sample_set, var):
    '''
    计算相邻评分的中位数，以便进行决策树二元切分
    param sample_set: 待切分样本
    param var: 分割变量名称
    '''
    var_list = list(np.unique(sample_set[var]))
    var_median_list = []
    for i in range(len(var_list) -1):
        var_median = (var_list[i] + var_list[i+1]) / 2
        var_median_list.append(var_median)
    return var_median_list

def choose_best_split(sample_set, var, min_sample):
    '''
    使用CART分类决策树选择最好的样本切分点
    返回切分点
    param sample_set: 待切分样本
    param var: 分割变量名称
    param min_sample: 待切分样本的最小样本量(限制条件)
    '''
    # 根据样本评分计算相邻不同分数的中间值
    target = "LS_type"
    score_median_list = calc_score_median(sample_set, var)
    median_len = len(score_median_list)
    #print("The number of unique values in feature",var,"is",median_len)
    sample_cnt = sample_set.shape[0]    # The number of sample 
    sample1_cnt = sum(sample_set[target])  # The number of class 1
    sample0_cnt =  sample_cnt- sample1_cnt   # The number of class 2
    Gini = 1 - np.square(sample1_cnt / sample_cnt) - np.square(sample0_cnt / sample_cnt)
    
    bestGini = 0.0; bestSplit_point = 0.0; bestSplit_position = 0.0
    
    for i in range(1,median_len,10000):
        #print("In the function choose_best_split, it is the",i,"split point")
        left = sample_set[sample_set[var] < score_median_list[i]]
        right = sample_set[sample_set[var] > score_median_list[i]]
        
        left_cnt = left.shape[0]; right_cnt = right.shape[0]
        left1_cnt = sum(left[target]); right1_cnt = sum(right[target])
        left0_cnt =  left_cnt - left1_cnt; right0_cnt =  right_cnt - right1_cnt
        left_ratio = left_cnt / sample_cnt; right_ratio = right_cnt / sample_cnt
        
        if left_cnt < min_sample or right_cnt < min_sample:
            continue
        
        Gini_left = 1 - np.square(left1_cnt / left_cnt) - np.square(left0_cnt / left_cnt)
        Gini_right = 1 - np.square(right1_cnt / right_cnt) - np.square(right0_cnt / right_cnt)
        Gini_temp = Gini - (left_ratio * Gini_left + right_ratio * Gini_right)
        if Gini_temp > bestGini:
            bestGini = Gini_temp; bestSplit_point = score_median_list[i]
            if median_len > 1:
                bestSplit_position = i / (median_len - 1)
            else:
                bestSplit_position = i / median_len
        else:
            continue
               
    Gini = Gini - bestGini
    return bestSplit_point, bestSplit_position

def bining_data_split(sample_set, var, min_sample, split_list):
    '''
    划分数据找到最优分割点list (split the data to find the optimal split point)
    param sample_set: 待切分样本 (dataset to be split)
    param var: 分割变量名称 (feature to be split)
    param min_sample: 待切分样本的最小样本量(限制条件) (the min volume of sample when spliting)
    param split_list: 最优分割点list  (The optimal split point)
    '''
    split, position = choose_best_split(sample_set, var, min_sample)
    if split != 0.0:
        split_list.append(split)
    # 根据分割点划分数据集，继续进行划分
    sample_set_left = sample_set[sample_set[var] < split]
    sample_set_right = sample_set[sample_set[var] > split]
    # 如果左子树样本量超过2倍最小样本量，且分割点不是第一个分割点，则切分左子树
    if len(sample_set_left) >= min_sample * 2 and position not in [0.0, 1.0]:
        #print("to left tree")
        bining_data_split(sample_set_left, var, min_sample, split_list)
    else:
        None
    # 如果右子树样本量超过2倍最小样本量，且分割点不是最后一个分割点，则切分右子树
    if len(sample_set_right) >= min_sample * 2 and position not in [0.0, 1.0]:
        #print("to the right")
        bining_data_split(sample_set_right, var, min_sample, split_list)
    else:
        None
        
def get_bestsplit_list(sample_set, var):
    '''
    根据分箱得到最优分割点list (Get the list of th optimal split points) 
    param sample_set: 待切分样本 (The dataset to be split)
    param var: 分割变量名称  (The name of feature to be split)
    '''
    # 计算最小样本阈值（终止条件）
    min_df = sample_set.shape[0] * 0.05
    split_list = []
    # 计算第一个和最后一个分割点
    bining_data_split(sample_set, var, min_df, split_list)
    return split_list

def ContinousToCategorical(data, features):
    d = data.copy()
    split_info = []
    all_splits = []
    for feature in features:
        print("Feature:",feature)
        feature_min = min(d[feature])-1   # min-1 and max+1 is in order to avoid missing data
        feature_max = max(d[feature])+1
        feature_split = get_bestsplit_list(d,feature)
        all_splits.append(feature_split)
        feature_split.append(feature_min)
        feature_split.append(feature_max)
        feature_split = np.sort(feature_split)
        a = pd.cut(d[feature], feature_split).cat.categories
        for i in range(len(a)):
            print("The",i+1,"class is in the interval",a[i])
        ind = a.astype(str)
        split_info.append(ind)
        d[feature] = pd.cut(d[feature], feature_split,labels=False)
        print("------------------------")
    return d, split_info, all_splits


def Outlier(data, variable="curvature", option=1): 
    if variable not in data:
        return
    data_copy = data.copy()
    if option == 1:
        variable_des = data[variable].describe()
        valid_max = variable_des["50%"] + 3*(variable_des['75%'] - variable_des["50%"])
        valid_min = variable_des["50%"] - 3*(variable_des["50%"] - variable_des["25%"])
        data_copy[variable] = data_copy[variable].clip(valid_min, valid_max)
        return data_copy
    else:
        data_WithoutOutliers = data_copy[np.abs(data_copy[variable]-data_copy[variable].mean())<=(3*data_copy[variable].std())] 
        return data_WithoutOutliers

def Data_processing_with_T_link(random_state):
    data_all = pd.read_csv('ROC_data.csv') 
    print("The size of original data is", len(data_all))
    data_use = data_all.loc[:,['LS_type','DA_area','slope_tan',
           'elevation', 'curvature', 'aspect', 'wet_index', 'litho', 'lulc']]

    data_use.loc[data_use.LS_type != 8, "LS_type" ] = 1
    data_use.loc[data_use.LS_type == 8, "LS_type" ] = 0

    # all areas
    data = data_use.drop(["DA_area"],axis=1)

    # log transformation to correct the skewness of data
    data["slope_tan"] = np.log(data["slope_tan"])

    #--------------------------Data cleaning----------------------------------------------------
    data = Outlier(data, "curvature", option=2)
    #data = Outlier(data,'slope_tan',option=2)
    #data = Outlier(data,"elevation",option=2)
    #data = Outlier(data,"aspect",option=2)
    #data = Outlier(data,"wet_index",option=2)
    print("The size of data after cleaning",len(data))
    print("There are %d data are deleted or be seen as outliers" % (len(data_all) - len(data)) )
    print("-------------------------------This part finished----------------------------------------------")
    #-------------------------------------------------------------------------------------------

    # split the whole data into 2 pieces: training and testing data
    y = data["LS_type"]
    y = np.ravel(y)
    X = data.drop(["LS_type"],axis=1)
    xtrain, xtest, ytrain, ytest = train_test_split(X, y, test_size=0.2, random_state=random_state)

    #--------------------------Tomek link-------------------------------------------------------
    tlink = TomekLinks(random_state=0,ratio="majority")
    X_resampled, y_resampled = tlink.fit_sample(xtrain, ytrain)
    X_resampled = pd.DataFrame(X_resampled) 
    col = data.columns[1:]
    X_resampled.columns = col
    print (sorted(Counter(y_resampled).items()))
    print("Tomek link algorithm has deleted %s sample from the majority class" % (len(ytrain)-len(y_resampled)))
    print("-------------------------------This part finished----------------------------------------------")
    #-------------------------------------------------------------------------------------------

    #---------------------Transform continuous features into categorical types-----------------
    cols = data.columns[1:6] # the columns needed to be transformed
    if sum(cols == ['slope_tan', 'elevation', 'curvature', 'aspect', 'wet_index']) != 5:
        print("mistankes: wrong columns!")
        return 
    
    # In practice, we do not need data_categorical but we need split_info, all_splits
    _ , split_info, all_splits = ContinousToCategorical(data, cols)

    # Training
    X_resampled["lulc"] = X_resampled["lulc"].apply(lambda x: ((x//10) -2))
    training = X_resampled.copy()
    training['slope_tan'] = pd.cut(training["slope_tan"], bins=sorted(all_splits[0]), labels=False, include_lowest=False)
    training['elevation'] = pd.cut(training["elevation"], bins=sorted(all_splits[1]), labels=False, include_lowest=False)
    training['curvature'] = pd.cut(training['curvature'], bins=sorted(all_splits[2]), labels=False, include_lowest=False)
    training['aspect'] = pd.cut(training["aspect"], bins=sorted(all_splits[3]), labels=False, include_lowest=False)
    training['wet_index'] = pd.cut(training["wet_index"], bins=sorted(all_splits[4]), labels=False, include_lowest=False)

    class_features = training.columns 
    training = pd.get_dummies(training, columns=class_features, drop_first=True)
    training.insert(loc=0, column='intercept', value=[1]*training.shape[0])

    # Testing
    xtest["lulc"] = xtest["lulc"].apply(lambda x: ((x//10) -2))
    testing = xtest.copy()
    testing['slope_tan'] = pd.cut(testing["slope_tan"], bins=sorted(all_splits[0]), labels=False, include_lowest=False)
    testing['elevation'] = pd.cut(testing["elevation"], bins=sorted(all_splits[1]), labels=False, include_lowest=False)
    testing['curvature'] = pd.cut(testing['curvature'], bins=sorted(all_splits[2]), labels=False, include_lowest=False)
    testing['aspect'] = pd.cut(testing["aspect"], bins=sorted(all_splits[3]), labels=False, include_lowest=False)
    testing['wet_index'] = pd.cut(testing["wet_index"], bins=sorted(all_splits[4]), labels=False, include_lowest=False)

    testing = pd.get_dummies(testing, columns=class_features, drop_first=True)
    testing.insert(loc=0, column='intercept', value=[1]*testing.shape[0])
    testing.index = range(len(testing))

    if training.shape[1] != testing.shape[1] or sum(training.columns == testing.columns) != testing.shape[1]:
        print("mistakes: the column index of training and testing are different")
        return

    print("-------------------------------This part finished----------------------------------------------")
    #-------------------------------------------------------------------------------------------
    return training,  y_resampled, testing, ytest


xtrain, ytrain, xtest, ytest =  Data_processing_with_T_link(random_state=951129)


# a check for training and testing
xtrain.head()
xtest.index = range(len(xtest))
xtest.head()

##############################################################




##############################################################
# STORE DATA IF NECESSARY
tt = xtrain.copy()
tt.insert(loc=0, column='LS_type', value=ytrain)
tt.to_csv('training.csv', header=True, index= False)
tt2 = xtest.copy()
tt2.insert(loc=0, column='LS_type', value=ytest)
tt2.to_csv('testing.csv', header=True, index= False)
trainings = pd.read_csv('training.csv')
testings = pd.read_csv('testing.csv')
xtrain = trainings.iloc[:,1:]
ytrain = trainings['LS_type']
xtest = testings.iloc[:, 1:]
ytest = testings['LS_type']
##############################################################



##############################################################

# TRAIN MODELS
# loss function: Weighted cross entropy

def get_roc(pos_prob,y_true):
    pos = y_true[y_true==1]
    neg = y_true[y_true==0]
    threshold = np.sort(pos_prob)[::-1]  # cutoff value include the range of all results
    y = y_true[pos_prob.argsort()[::-1]]
    tpr_all = [0] ; fpr_all = [0]
    tpr = 0 ; fpr = 0
    x_step = 1/float(len(neg))
    y_step = 1/float(len(pos))
    y_sum = 0                             
    for i in range(len(threshold)):
        if y[i] == 1:
            tpr += y_step
            tpr_all.append(tpr)
            fpr_all.append(fpr)
        else:
            fpr += x_step
            fpr_all.append(fpr)
            tpr_all.append(tpr)
            y_sum += tpr
    return tpr_all,fpr_all,y_sum*x_step

def AdjustClassWeight(p, training, ytrain, testing, ytest):
    prec = []
    r = []
    f1 = []
    acc = []
    AUC = []
    for weight in p:
        class_weight = {0:weight, 1:1-weight}
        #print("weight=",class_weight)
        model = LogisticRegression(class_weight=class_weight)
        model.fit(training, ytrain)
        prediction = model.predict(testing)
        accuracy = model.score(testing, ytest)
        acc.append(accuracy)
        f1score = f1_score(ytest, prediction, average='binary')
        f1.append(f1score)
        precision = precision_score(ytest, prediction, average='binary')
        prec.append(precision) 
        recall = recall_score(ytest, prediction , average='binary')
        r.append(recall)
        prob_prediction = model.predict_proba(testing)[:,1]
        _,_,auc = get_roc(prob_prediction,ytest)  
        AUC.append(auc)
        d = pd.DataFrame({"Precision":prec,"Recall":r, "F1 Score":f1,"AUC":AUC, "Accuracy":acc})
    return d

# Plot and optimized
p = np.arange(0.01,0.35,0.01)
d = AdjustClassWeight(p, training, ytrain, testing, ytest)  
# plot f1 score and auc for different weights
d["weight of '0'"] = p
d = d[:30]
plt.plot(d["weight of \'0\'"], d["F1 Score"],linestyle='dashed', color='k', marker='o' )
plt.xlabel("weight of class 0")
plt.ylabel("F1 score value")
plt.title("F1 Score under different weights")
plt.show()
plt.plot(d["weight of \'0\'"], d["Precision"],linestyle='dashed', color='blue', marker='*', label="Precision")
plt.plot(d["weight of \'0\'"], d["Recall"],linestyle='dashed', color='red', marker='o', label="Recall" )
plt.legend()
plt.xlabel("weight of class 0")
plt.title("Precision and Recall under different weights")
ind_opt = d["F1 Score"].idxmax()
print("the optimal weight is",{0:p[ind_opt],1:1-p[ind_opt]})
d.loc[ind_opt]

model = LogisticRegression(class_weight={0:0.08,1:0.92})
model.fit(xtrain, ytrain)
y_predprob = model.predict_proba(xtrain)[:,1]

# Training
print("------------------Training performance--------------------------------------")
print ("AUC Score (test): %f" % roc_auc_score(ytrain, y_predprob))
prediction =  model.predict(xtrain)
print("Precison is",precision_score(ytrain, prediction, average='binary'))
print("Recall is",recall_score(ytrain, prediction , average='binary'))
print("F1 score is", f1_score(ytrain, prediction , average='binary'))
print("Accuracy is",  model.score(xtrain, ytrain))
print(confusion_matrix(ytrain, prediction))
print("-------------------------Finished--------------------------------------------")


# Testing
print("------------------Testing performance--------------------------------------")
y_predprob = model.predict_proba(xtest)[:,1]
print ("AUC Score (test): %f" % roc_auc_score(ytest, y_predprob))
prediction =  model.predict(xtest)
print("Precison is",precision_score(ytest, prediction, average='binary'))
print("Recall is",recall_score(ytest, prediction , average='binary'))
print("F1 score is", f1_score(ytest, prediction , average='binary'))
print("Accuracy is",  model.score(xtest, ytest))
print(confusion_matrix(ytest, prediction))
print("-------------------------Finished--------------------------------------------")


##############################################################
# FOCAL LOSS
from keras.models import Sequential
import tensorflow as tf
#from tensorflow import keras
from keras.layers import Dense, Dropout, Activation, Flatten
from keras import backend as K
from keras.callbacks import Callback

def f1(y_true, y_pred):
    def recall(y_true, y_pred):
        """Recall metric.
        Only computes a batch-wise average of recall.
        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall
 
    def precision(y_true, y_pred):
        """Precision metric.
        Only computes a batch-wise average of precision.
        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

def binary_focal_loss(gamma=2.0, alpha=0.25):
    """
    Implementation of Focal Loss from the paper in multiclass classification
    Formula:
        loss = -alpha_t*((1-p_t)^gamma)*log(p_t)
        
        p_t = y_pred, if y_true = 1
        p_t = 1-y_pred, otherwise
        
        alpha_t = alpha, if y_true=1
        alpha_t = 1-alpha, otherwise
        
        cross_entropy = -log(p_t)
    Parameters:
        alpha -- the same as wighting factor in balanced cross entropy
        gamma -- focusing parameter for modulating factor (1-p)
    Default value:
        gamma -- 2.0 as mentioned in the paper
        alpha -- 0.25 as mentioned in the paper
    """
    def focal_loss(y_true, y_pred):
        # Define epsilon so that the backpropagation will not result in NaN
        # for 0 divisor case
        
        #epsilon = K.epsilon()
        
        # Add the epsilon to prediction value
        #y_pred = y_pred + epsilon
        # Clip the prediciton value
        
        #y_pred = K.clip(y_pred, epsilon, 1.0-epsilon)
        
        # Calculate p_t
        p_t = tf.where(K.equal(y_true, 1), y_pred, 1-y_pred)
        # Calculate alpha_t
        alpha_factor = K.ones_like(y_true)*alpha
        alpha_t = tf.where(K.equal(y_true, 1), alpha_factor, 1-alpha_factor)
        # Calculate cross entropy
        cross_entropy = -K.log(p_t)
        weight = alpha_t * K.pow((1-p_t), gamma)
        # Calculate focal loss
        loss = weight * cross_entropy
        # Sum the losses in mini_batch
        loss = K.sum(loss, axis=1)
        
        
        return loss*1000
    
    return focal_loss 

model = Sequential()

input_dim = xtrain.shape[1]
#nb_classes = y_train.shape[1]

model = Sequential()
model.add(Dense(input_dim=input_dim, units=1))
model.add(Activation('sigmoid'))
from keras import optimizers
opt = optimizers.Nadam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=1e-08, schedule_decay=0.004)
model.compile(loss=binary_focal_loss(gamma=2,alpha=0.9),
              optimizer=opt,metrics=[f1])

history = model.fit(xtrain, ytrain, epochs=15, batch_size=500, verbose=1) 


# Training
print("------------------Training performance--------------------------------------")
print ("AUC Score (test): %f" % roc_auc_score(ytrain, y_predprob))
prediction =  model.predict(xtrain)
print("Precison is",precision_score(ytrain, prediction, average='binary'))
print("Recall is",recall_score(ytrain, prediction , average='binary'))
print("F1 score is", f1_score(ytrain, prediction , average='binary'))
print("Accuracy is",  model.score(xtrain, ytrain))
print(confusion_matrix(ytrain, prediction))
print("-------------------------Finished--------------------------------------------")


# Testing
print("------------------Testing performance--------------------------------------")
y_predprob = model.predict_proba(xtest)[:,1]
print ("AUC Score (test): %f" % roc_auc_score(ytest, y_predprob))
prediction =  model.predict(xtest)
print("Precison is",precision_score(ytest, prediction, average='binary'))
print("Recall is",recall_score(ytest, prediction , average='binary'))
print("F1 score is", f1_score(ytest, prediction , average='binary'))
print("Accuracy is",  model.score(xtest, ytest))
print(confusion_matrix(ytest, prediction))
print("-------------------------Finished--------------------------------------------")










