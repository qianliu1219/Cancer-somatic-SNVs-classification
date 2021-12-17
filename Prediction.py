# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:40:40 2019

@author: Nikta
"""
#Predicting pathogenicity of 1) Metabrik(coding) 2)TCGA (coding and non_coding)
# two models: 1)svm-rfb-normal data coding and non-coding 2)logreg imput data coding and non-coding
import pandas as pd
from sklearn.utils import shuffle# Mixing the pos and neg data by shuffeling:
from statsmodels.imputation.mice import MICEData # Using this pakage for imputing cnumeric features
from sklearn import preprocessing
#from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC 
import numpy as np
import matplotlib.pyplot as plt
         
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
    
coding=readdata("fin_af_cmb_cd2nd.csv","fin_af_neg_cd2nd.csv")
coding=coding[~coding.index.duplicated(keep=False)]
#coding=coding.iloc[0:1000,:]#test
non_coding=readdata("fin_af_cmb_noncd2nd.csv","fin_af_neg_noncd2nd.csv")
non_coding=non_coding[~non_coding.index.duplicated(keep=False)]
#non_coding=non_coding.iloc[0:2000,:]#test

data_names=[coding,non_coding]
imp=[]
norm=[]
for i in data_names:
    imp_data=imputation(i)
    imp.append(imp_data)
    norm.append(normalization(imp_data))
    
data={"cd":[imp[0],norm[0]],"ncd":[imp[1],norm[1]]}

#Prediction data:
Metabrik=pd.read_csv("metabrik2nd.csv",index_col=0)# duplication is checked before! No duplications!
Metabrik_imp=imputation(Metabrik)
Metabrik_norm=normalization(Metabrik_imp)
#Metabrik.head()
#Metabrik.isnull().sum()
#Metabrik=Metabrik.dropna(axis=0)#0, or ‘index’ : Drop rows which contain missing values.

TCGA_cd=pd.read_csv("TCGA_CODING2nd.csv",index_col=0)
TCGA_cd_imp=imputation(TCGA_cd)
TCGA_cd_norm=normalization(TCGA_cd_imp)

TCGA_Ncd=pd.read_csv("TCGA_NonCODING2nd.csv",index_col=0)
TCGA_Ncd_imp=imputation(TCGA_Ncd)
TCGA_Ncd_norm=normalization(TCGA_Ncd_imp)

#Models:
Yc=data["cd"][0].Lable
Ync=data["ncd"][0].Lable

Xc_imp=data["cd"][0].iloc[:,0:-1]#for lasso
Xc_norm=data["cd"][1].iloc[:,0:-1]#for svm
Xnc_imp=data["ncd"][0].iloc[:,0:-1]#for lasso
Xnc_norm=data["ncd"][1].iloc[:,0:-1]#for svm

#SVM Models:
clf= SVC(kernel='rbf',gamma='auto', probability=1) #Coding data:(predcitng metbrik + coding TCGA)
SVM_cd=clf.fit(Xc_norm, Yc)# Training for normalized coding data

svm_met_prob=SVM_cd.predict_proba(Metabrik_norm)[:, 1]# This would show the probabilities
svm_met_lab=SVM_cd.predict(Metabrik_norm) # This would show the prdicted lables

svm_cd_tcga_prob=SVM_cd.predict_proba(TCGA_cd_norm)[:, 1]
svm_cd_tcga_lab=SVM_cd.predict(TCGA_cd_norm)

SVM_Ncd=clf.fit(Xnc_norm , Ync)# Training for normalized non-coding data and predicting TCGA Noncoding

svm_Nc_tcga_prob=SVM_Ncd.predict_proba(TCGA_Ncd_norm)[:, 1]
svm_Nc_tcga_lab=SVM_Ncd.predict(TCGA_Ncd_norm)

#Logestic Models:
#logreg = LogisticRegression(penalty="l1", solver='liblinear')
#lasso_cd=logreg.fit(Xc_imp,Yc)

#lasso_met_prob=lasso_cd.predict_proba(Metabrik_imp)[:, 1]# This would show the probabilities
#lasso_met_lab=lasso_cd.predict(Metabrik_imp) # This would show the prdicted lables

#lasso_cd_tcga_prob=lasso_cd.predict_proba(TCGA_cd_imp)[:, 1]
#lasso_cd_tcga_lab=lasso_cd.predict(TCGA_cd_imp)

#lasso_Ncd=logreg.fit(Xnc_imp,Ync)
#lasso_Nc_tcga_prob=lasso_Ncd.predict_proba(TCGA_Ncd_imp)[:, 1]
#lasso_Nc_tcga_lab=lasso_Ncd.predict(TCGA_Ncd_imp)

################################################################### 
preicted_res_prob=[svm_met_prob,svm_cd_tcga_prob,svm_Nc_tcga_prob]#,lasso_met_prob,lasso_cd_tcga_prob,lasso_Nc_tcga_prob]
#predivted_lables_default=[svm_met_lab,svm_cd_tcga_lab,svm_Nc_tcga_lab]#,lasso_met_lab,lasso_cd_tcga_lab,lasso_Nc_tcga_lab]#threshold is 0.5!
#thresholds=[0.55,0.55,0.41,0.43,0.43,0.34]#choose based on each data-set and model respectively!
Metabrik_id_gene=pd.read_csv("METABRIK_Genes.csv", index_col=0)
TCGA_CD_id_gene=pd.read_csv("TCGA_CODING_Genes.csv", index_col=0)
TCGA_NCD_id_gene=pd.read_csv("TCGA_NonCODING_Genes.csv", index_col=0)

org_metabrik_45=pd.read_csv("Metabrik_mutation_young.csv") # young metabrik patients, with 38 assembly
org_metabrik_45["Mut_ID"]= org_metabrik_45["Chromosome"].astype(str) + ":"+ org_metabrik_45["38"].astype(str) + "_"+ org_metabrik_45["Reference_Allele"] + ">"+ org_metabrik_45["Tumor_Seq_Allele2"]# this is id based on 38 assembly 
org_tcga_mixed=pd.read_csv("tcga_young.csv")# young tcga data
org_tcga_mixed["Mut_ID"]= org_tcga_mixed["Chromosome"].astype(str) + ":"+ org_tcga_mixed["pos38"].astype(str) + "_"+ org_tcga_mixed["Reference_Allele"] + ">"+ org_tcga_mixed["Tumor_Seq_Allele2"]# this is id based on 38 assembly 

def gene_snp(gene_data,pred_data):# index of gene-data should be the snp ids!
    dic={}
    for k in  gene_data.index:      
        if k in pred_data:
           dic_key=gene_data.loc[k,"GeneName"]
           if dic_key not in dic.keys():
              dic[dic_key]=[]
           dic[dic_key].append(k)  
    return(dic)          


def freq(org_data,pred_data,t=0): # Counting the frequency of pos predictions accros whole data set, t is threshold
    if t==0:
       freq_dic={} 
       for ID in pred_data:
           freq_dic[ID]=0
           for idd in org_data:
               if ID==idd:
                   freq_dic[ID]=freq_dic[ID]+1
       return(freq_dic)    # showing each id along with its frequncy              
    else:  
        count_list=[] # count of the number of snps bigger than threshold
        for ID in pred_data:
            c=0
            for idd in org_data:      
                if ID==idd:
                   c=c+1  
            if c>=t:
               count_list.append(c)  
        return(len(count_list))   # returning the number of ids with the frequency of bigger than the threshold! 

def gene_finder(gene_data,org_data,pred_data,t):
   dic={}
   for ID in pred_data:
            c=0
            for idd in org_data:      
                if ID==idd:
                   c=c+1  
            if c>=t:
               dic_key=gene_data.loc[ID,"GeneName"]
               if dic_key not in dic.keys():
                  dic[dic_key]=[]
               dic[dic_key].append(ID) 
   return(dic) 


thresholds=[[0.55,0.55,0.41],[0.6,0.6,0.45],[0.65,0.65,0.5],[0.7,0.7,0.55],[0.75,0.75,0.6],[0.8,0.8,0.65],[0.85,0.85,0.7]]
count_positive_predictions=[[],[],[],[],[],[],[]]
#positive_lables_defualt=[]
index_positive_pedictions=[[],[],[],[],[],[],[]]
IDs_positive_predictions=[[],[],[],[],[],[],[]]
gene_snp_dics=[[],[],[],[],[],[],[]]
count_affected_genes= [[],[],[],[],[],[],[]]
snp_reccurency_dics=[[],[],[],[],[],[],[]]
count_recurrency_2=[[],[],[],[],[],[],[]]
count_recurrency_3=[[],[],[],[],[],[],[]]
count_recurrency_4=[[],[],[],[],[],[],[]]
genes_rec_2=[[],[],[],[],[],[],[]]
genes_rec_3=[[],[],[],[],[],[],[]]
genes_rec_4=[[],[],[],[],[],[],[]]
count_genes_rec_2=[[],[],[],[],[],[],[]]
count_genes_rec_3=[[],[],[],[],[],[],[]]
count_genes_rec_4=[[],[],[],[],[],[],[]]

for ii in range(len(thresholds)):
    
    for i in range(len(preicted_res_prob)):
        thresh=thresholds[ii][i]
        count_positive_predictions[ii].append(np.sum(preicted_res_prob[i] > thresh))
        #positive_lables_defualt.append(np.sum(predivted_lables_default[i]))
        index_positive_pedictions[ii].append([j for j in np.where(preicted_res_prob[i]> thresh)])
        if i==0:
           data=Metabrik
           whole_data=Metabrik_id_gene
           original_data=org_metabrik_45
        elif i==1:
           data= TCGA_cd_imp
           whole_data=TCGA_CD_id_gene
           original_data=org_tcga_mixed
        else:
           data= TCGA_Ncd_imp
           whole_data= TCGA_NCD_id_gene
           original_data=org_tcga_mixed
           
        prediction_data=data.iloc[index_positive_pedictions[ii][i][0]].index.tolist() #for up_coming functions! 
        IDs_positive_predictions[ii].append(prediction_data)# it works with booliean array!^-^ for svm
        
        gene_snp_dic=gene_snp(whole_data,prediction_data)# reporting the genes and the snps they have been affected by
        gene_snp_dics[ii].append(gene_snp_dic)
        count_affected_genes[ii].append(len(gene_snp_dic)) # reporting the total number of genes affeced by positive snps! 
       
        snp_reccurency_dic=freq(original_data["Mut_ID"],prediction_data) # frequency of each pos prediction accross originL young DATA      
        snp_reccurency_dics[ii].append(snp_reccurency_dic)
        count_recurrency_2[ii].append(freq(original_data["Mut_ID"],prediction_data,2))
        count_recurrency_3[ii].append(freq(original_data["Mut_ID"],prediction_data,3))
        count_recurrency_4[ii].append(freq(original_data["Mut_ID"],prediction_data,4))
        genes_rec_2[ii].append(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,2))
        genes_rec_3[ii].append(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,3))
        genes_rec_4[ii].append(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,4))
        count_genes_rec_2[ii].append(len(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,2)))
        count_genes_rec_3[ii].append(len(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,3)))
        count_genes_rec_4[ii].append(len(gene_finder(whole_data,original_data["Mut_ID"],prediction_data,4)))


#selina_genes=pd.read_csv("selina_gense.csv")
# Data frames containing genes along with the SNVs affected
Gene_SNV_met=pd.DataFrame.from_dict(gene_snp_dics[0][0], orient='index')
Gene_SNV_met.to_csv("Gene_SNVPosPred_met.csv")# showing the snvs each gene is affcted by

Gene_SNV_TCGA_cd=pd.DataFrame.from_dict(gene_snp_dics[0][1], orient='index')
Gene_SNV_TCGA_cd.to_csv("Gene_SNVPosPred_TCGA_cd.csv")# showing the snvs each gene is affcted by

Gene_SNV_TCGA_Ncd=pd.DataFrame.from_dict(gene_snp_dics[0][2], orient='index')
Gene_SNV_TCGA_Ncd.to_csv("Gene_SNVPosPred_TCGA_Ncd.csv")# showing the snvs each gene is affcted by

pd.DataFrame(IDs_positive_predictions[5][0]).to_csv("Met_500_snv.csv")
pd.DataFrame(IDs_positive_predictions[5][1]).to_csv("TCGA_cd_1860.csv")
pd.DataFrame(IDs_positive_predictions[5][2]).to_csv("TCGA_ncd_435.csv")


#plotting the results
x_c=[0.55,0.6,0.65,0.7,0.75,0.8,0.85]
snv_met=[]
gene_met=[]
for i in range(7):
    snv_met.append(count_positive_predictions[i][0])
    gene_met.append(count_affected_genes[i][0])

#results from server:
snv_met=[960, 880, 799, 705, 608, 502, 369]
gene_met=[155, 148, 144, 138, 130, 117, 105]
    

fig,ax1=plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Thresholds')
ax1.set_ylabel('Count of positive predicted SNVs', color=color)
ax1.plot(x_c, snv_met, 'ro', x_c, snv_met,'k')
ax1.tick_params(axis='y', labelcolor=color)
ax2=ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel("count of affected genes", color=color)
ax2.plot(x_c, gene_met, 'bo', x_c, gene_met, 'k')
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
plt.savefig("gene_snv_freq_met")

snv_TCGA_CD=[]
gene_TCGA_CD=[]
for i in range(7):
    snv_TCGA_CD.append(count_positive_predictions[i][1])
    gene_TCGA_CD.append(count_affected_genes[i][1])
#results from server:
snv_TCGA_CD=[2183, 1988, 1807, 1608, 1390, 1110, 829]
gene_TCGA_CD=[1651, 1519, 1395, 1260, 1104, 907, 703]    
    
fig,ax1=plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Thresholds')
ax1.set_ylabel('Count of positive predicted SNVs', color=color)
ax1.plot(x_c, snv_TCGA_CD, 'ro', x_c, snv_TCGA_CD,'k')
ax1.tick_params(axis='y', labelcolor=color)
ax2=ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel("count of affected genes", color=color)
ax2.plot(x_c, gene_TCGA_CD, 'bo', x_c, gene_TCGA_CD, 'k')
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
plt.savefig("gene_snv_freq_TCGA_CD")    

x_nc=[0.41,0.45,0.5,0.55,0.6,0.65,0.7]
snv_TCGA_NCD=[] 
gene_TCGA_NCD=[]

for i in range(7):
    snv_TCGA_NCD.append(count_positive_predictions[i][2])
    gene_TCGA_NCD.append(count_affected_genes[i][2])

#results from server:
snv_TCGA_NCD=[476, 431, 366, 312, 263, 206, 164]
gene_TCGA_NCD=[174, 165, 150, 141, 135, 123, 112]    
     
fig,ax1=plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Thresholds')
ax1.set_ylabel('Count of positive predicted SNVs', color=color)
ax1.plot(x_nc, snv_TCGA_NCD, 'ro', x_nc, snv_TCGA_NCD,'k')
ax1.tick_params(axis='y', labelcolor=color)
ax2=ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel("count of affected genes", color=color)
ax2.plot(x_nc, gene_TCGA_NCD, 'bo', x_nc, gene_TCGA_NCD, 'k')
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show() 
plt.savefig("gene_snv_freq_TCGA_NCD")
#########################################recurrence plots:
rec=[1,2,3,4,5,6,7,8,9]
rec_snv_met=[]
for jj in rec:
    rec_snv_met.append(freq(org_metabrik_45["Mut_ID"],IDs_positive_predictions[0][0],t=jj))# for the threshold of 0.55

























