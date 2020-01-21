# Self-defined Focal Loss functions

def Focal_Logistic_MiniBatch(X,y,alpha=1,gamma=0,class_weight={0:0.5,1:0.5},random_state=1):
    '''
    Gradient Descending: slow
    alpha: the learning rate
    X: training data
    y: training labels
    gamma: focusing parameter in focal loss, and if gamma=0, focal loss is the cross entropy
    class_weight: the weight of two different classes
    '''
    X = np.array(X)
    np.random.seed(random_state)
    alpha = alpha
    beta = np.random.randn(X.shape[1])
    weight_0 = class_weight[0]  # weight_0 is the class weight of y=0
    weight_1 = class_weight[1]  # weight_1 is the class weight of y=1 
    losses = []
    auc =[]
    batch_size = 1
    
    for T in range(100):
        prob = np.array(1./ (1 + np.exp(-np.matmul(X, beta)))).ravel()   # calculate P[Y=1|X,beta]
        prob_y = list(zip(prob, y))
        FLoss = -sum([weight_1*((1-p)**gamma)*np.log(p) if y == 1 else weight_0*(p**gamma)*np.log(1 - p) for p, y in prob_y]) / len(y)
        losses.append(FLoss)
        #auc_val = roc_auc_score(y,prob)
        #auc.append(auc_val)
        
        if T % 2 == 0:
            #print("T=" + str(T) + "Focal Loss=" + str(FLoss) + "AUC=" + str(auc_val))
            print("T=" + str(T) + "Focal Loss=" + str(FLoss))
            
        # calculate the derivatives
        deriv = np.zeros(X.shape[1])
        for i in range(batch_size):
            
            idxs = np.random.randint(0,X.shape[0], size=batch_size)
            X_temp = X.take(idxs, axis=0)
            y_temp = y.take(idxs, axis=0)
            
            neg_exp = np.array(np.exp(-np.matmul(X_temp[i,:], beta))).ravel()
            
            pos_exp = np.array( np.exp(+np.matmul(X_temp[i,:], beta))).ravel() 
            
            kernel_val_1 = weight_0 * (1-y_temp[i]) * (1. / (1 + neg_exp)**(gamma+1)) * (gamma*neg_exp*np.log(1+pos_exp) + 1)
            
            kernel_val_2 = weight_1 * y_temp[i] * (1. / (1 + pos_exp)**(gamma+1)) * (gamma*pos_exp*np.log(1+neg_exp) + 1)
            
            deriv += np.asarray(X_temp[i,:]).ravel() * (kernel_val_1 - kernel_val_2) 
            
        
        deriv /= len(y)
        beta -= alpha * deriv 
    
    return beta

def Predict_Proba_Focal_Logistic(X, beta):
    prob = np.array(1. / (1 + np.exp(-np.matmul(X, beta)))).ravel()
    return prob
    
def Predict_label_Focal_Logistic(X, beta, cutoff_val=0.5):
    prob = np.array(1. / (1 + np.exp(-np.matmul(X, beta)))).ravel()
    prediction =  1*(prob>cutoff_val)
    return prediction