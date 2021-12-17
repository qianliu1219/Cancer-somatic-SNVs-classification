# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:26:06 2019

@author: SONY
"""

#collecting negative samples: (label as negative in the final model)
import pandas as pd
import numpy as np
Neg=pd.read_csv("CosmicMutantExport.tsv",error_bad_lines=False, sep="\t")
Neg_arr=np.array(Neg)
c=0
for i in range(len(Neg_arr)):
     if Neg_arr[i,22]==38.0 and ">" in Neg_arr[i,17] and Neg_arr[i,25] == "y" and Neg_arr[i,29]=="Confirmed somatic variant" and "Substitution" in  Neg_arr[i,19] : 
         c=c+1
fin_Neg=np.array(np.random.rand (c,35),dtype='str')
j=0
for i in range(len(Neg_arr)):
     if  Neg_arr[i,22]==38.0 and ">" in Neg_arr[i,17] and Neg_arr[i,25] =="y" and Neg_arr[i,29]=="Confirmed somatic variant" and "Substitution" in Neg_arr[i,19]   :
         fin_Neg[j,:]=Neg_arr[i,:]
         j=j+1
fin_Neg_df=pd.DataFrame (fin_Neg, columns=Neg.columns)
fin_Neg_df=fin_Neg_df.drop_duplicates(subset=['Mutation CDS','Mutation genome position','Primary site'] , keep = 'first' , inplace=False)
fin_Neg_df["ID_for_match_wdbsnp"]=fin_Neg_df['Mutation genome position']+ "_" + fin_Neg_df['Mutation CDS'].str[-3] + ">" + fin_Neg_df['Mutation CDS'].str[-1] 


# By now we have the negative data from cosmic. 
#Then I will open the data from dbsnp and try to match these two sets of data.
#dbsnp= pd.read_csv("common_all_20180418.vcf", error_bad_lines=False, sep="\t")

dbsnp= open("common_all_20180418.vcf", "r")
for i in range(56):
       print(dbsnp.readline())
head= dbsnp.readline().strip().split("\t")
line = dbsnp.readline().strip().split("\t")
w=open("dbsnp_commonvar.csv","w")
w.write(",".join(map(str,head))+'\n')
i=0
while len(line)>1:
     if "COMMON=1" in line[-1].strip().split(";") and len(line[3])==1 and len(line[4])==1 and line[0] != "X" and line[0] != "Y" :
           newline=",".join(map(str,line[0:6]))+'\n'
           w.write(newline)
     line= dbsnp.readline().strip().split("\t")
     i=i+1
     if i%100000==0:
          print(i)  
dbsnp.close
w.close
final_dbsnp= pd.read_csv("dbsnp_commonvar.csv", error_bad_lines=False)
final_dbsnp["dbsnp_position_ID"]=final_dbsnp['#CHROM'].astype(str)+ ":" +  final_dbsnp['POS'].astype(str)+"-"+final_dbsnp['POS'].astype(str)+ "_" + final_dbsnp['REF']+">"+ final_dbsnp['ALT']
match=list(fin_Neg_df["ID_for_match_wdbsnp"].isin(final_dbsnp["dbsnp_position_ID"]))
indic=[i for i in range(len(match)) if match[i]]
Neg_cd_fin= fin_Neg_df.iloc[indic,:]
Neg_cd_fin.to_csv("Neg_cd_cosdbsnp.csv", sep="\t")
