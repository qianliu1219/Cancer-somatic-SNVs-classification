# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:53:45 2019

@author: Nikta
"""

# Modles funcions:
import pandas as pd
from sklearn.utils import shuffle# Mixing the pos and neg data by shuffeling:
from statsmodels.imputation.mice import MICEData # Using this pakage for imputing cnumeric features
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC 
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt 
from sklearn.metrics import roc_auc_score


def readdata(pos,neg):
    positive= pd.read_csv(pos,index_col=0)
    Negative=pd.read_csv(neg,index_col=0)
    data=pd.concat([positive, Negative], axis=0,ignore_index=0)
    data= shuffle(data)
    #print(data.isnull().sum())
    return(data)
    
def imputation(dada):
    imp=MICEData(dada)
    dada_imp=imp.data
    dada_imp.set_index(dada.index,inplace=True)
    return(dada_imp)

def normalization(data):
    data_norm=pd.DataFrame(preprocessing.scale(data))
    data_norm.columns=data.columns
    data_norm.set_index(data.index,inplace=True)
    return(data_norm)

def logistic(y,x,crossvalidation,penal='l1',cross=10):
    logreg = LogisticRegression(penalty=penal, solver='liblinear')#penalty='l1' means we are using lasso!
    if crossvalidation=="yes":
       scoring = {'acc': 'accuracy' , 'AUC': 'roc_auc'} 
       k_fold = KFold(n_splits=cross,random_state=50, shuffle=False)
       scores = cross_validate(logreg, x, y , scoring=scoring, cv=k_fold, return_train_score=True)
       return(scores['train_AUC'].mean(),scores['train_acc'].mean(),scores['test_AUC'].mean(), scores['test_acc'].mean())
    else:   
       X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3,random_state=50)
       logreg.fit(X_train,y_train)
       probs = logreg.predict_proba(X_test)[:, 1] 
       auc_onethird = roc_auc_score(y_test, probs)
       accuracy=logreg.score(X_test, y_test)
       return(auc_onethird,accuracy)
  
def svm(y,x,ker,crossvalidation,cross=10):
       if crossvalidation=="yes":
          #print("Im using the first loop")
          clf= SVC(kernel=ker,verbose=True)
          scoring = {'acc': 'accuracy' , 'AUC': 'roc_auc'}          
          k_fold = KFold(n_splits=cross,random_state=50, shuffle=True)
          scores = cross_validate(clf, x , y , scoring=scoring, cv=k_fold, return_train_score=True)
          return(scores['train_AUC'].mean(),scores['train_acc'].mean(),scores['test_AUC'].mean(), scores['test_acc'].mean())   
       else:
           X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3,random_state=50)
           #print("first line is done")
           clf = SVC(kernel=ker, probability=1).fit(X_train, y_train)
           #print("im done training")
           probs = clf.predict_proba(X_test)[:, 1] 
           auc_onethird = roc_auc_score(y_test, probs)
           accuracy=clf.score(X_test, y_test)
           return(auc_onethird,accuracy)
                     
    

coding=readdata("fin_af_cmb_cd2nd.csv","fin_af_neg_cd2nd.csv")
coding=coding.iloc[0:1000,:]#for test
coding=coding[~coding.index.duplicated(keep=False)]#new
non_coding=readdata("fin_af_cmb_noncd2nd.csv","fin_af_neg_noncd2nd.csv")
non_coding=non_coding.iloc[0:2000,:]#for test
non_coding=non_coding[~non_coding.index.duplicated(keep=False)]#new
coding_drop=coding.drop(["SIFTval","PolyPhenVal","Grantham","GerpRS", "GerpRSpval"],axis=1)
non_coding_drop=non_coding.drop(["GerpN","GerpS","EncodetotalRNA-max"],axis=1)
data_names=[coding,non_coding,coding_drop,non_coding_drop]
imp_files=[]
norm_files=[]
#Files for logreg model:
imp_cv_reg={}# 0 is coding, 1 is non_coding, 2 is coding_drop, 3 is non_coding drop
norm_cv_reg={}# 0 is coding, 1 is non_coding,2 is coding_drop, 3 is non_coding drop
imp_third_reg={}# 0 is coding, 1 is non_coding,2 is coding_drop, 3 is non_coding drop
norm_third_reg={}# 0 is coding, 1 is non_coding,2 is coding_drop, 3 is non_coding drop
#Files for svm model:
kernel=["rbf","sigmoid"]
imp_cv_svm=[{},{},{},{}]#first dic is coding and the second is non_coding
imp_third_svm=[{},{},{},{}]
norm_cv_svm=[{},{},{},{}]#first dic is coding and the second is non_coding
norm_third_svm=[{},{},{},{}]#first dic is coding and the second is non_coding
    
for i in data_names:
    imp_data=imputation(i)
    imp_files.append(imp_data)
    norm_files.append(normalization(imp_data))
   
for i in range(len(data_names)):#logreg model
    imp_cv_reg[i]=logistic(imp_files[i].Lable,imp_files[i].iloc[:,0:-1],crossvalidation="yes")#imputation and cv
    imp_third_reg[i]=logistic(imp_files[i].Lable,imp_files[i].iloc[:,0:-1],crossvalidation="no")#imputation without cv
    norm_cv_reg[i]=logistic(imp_files[i].Lable,norm_files[i].iloc[:,0:-1],crossvalidation="yes")#normalization and cv
    norm_third_reg[i]=logistic(imp_files[i].Lable,norm_files[i].iloc[:,0:-1],crossvalidation="no")#normalization without cv

for k in range(len(imp_cv_svm)):  
   print(k)    
   for ii in range(len(data_names)):
        print(ii)
        for j in kernel:
            print(j)
            imp_cv_svm[k][j]=svm(imp_files[ii].Lable,imp_files[ii].iloc[:,0:-1],j,crossvalidation="yes",cross=10)#imputation and cv
            imp_third_svm[k][j]=svm(imp_files[ii].Lable,imp_files[ii].iloc[:,0:-1],j,crossvalidation="no")##imputation without cv
            norm_cv_svm[k][j]=svm(imp_files[ii].Lable,norm_files[ii].iloc[:,0:-1],j,crossvalidation="yes",cross=10,)#normaization and cv
            norm_third_svm[k][j]=svm(imp_files[ii].Lable,norm_files[ii].iloc[:,0:-1],j,crossvalidation="no")#normaization without cv

   
    
    
    
 

           
k=0
j="rbf"
imp_cv_svm=[]            
kernel = ["rbf"]
for k in range(1):  
   print(k)    
   for ii in range(1):
        print(ii)
        for j in kernel:
            print(j)
            imp_cv_svm[k][j]=svm(imp_files[ii].Lable,imp_files[ii].iloc[:,0:-1],j,crossvalidation="yes",cross=10)#imputation and cv
            print("step 1 done")
            imp_third_svm[k][j]=svm(imp_files[ii].Lable,imp_files[ii].iloc[:,0:-1],j,crossvalidation="no")##imputation without cv
            print("step 2 done")
            norm_cv_svm[k][j]=svm(imp_files[ii].Lable,norm_files[ii].iloc[:,0:-1],j,crossvalidation="yes",cross=10,)#normaization and cv
            print("step 3 done")
            norm_third_svm[k][j]=svm(imp_files[ii].Lable,norm_files[ii].iloc[:,0:-1],j,crossvalidation="no")#normaization without cv
            print("step 4 done")

imp_cv_svm[0]["rbf"]=svm(imp_files[0].Lable,imp_files[0].iloc[:,0:-1],"rbf","yes")
imp_third_svm[0]["rbf"]=svm(imp_files[0].Lable,imp_files[0].iloc[:,0:-1],"rfb",crossvalidation="no")##imputation without cv
            



















































































     
            
            

