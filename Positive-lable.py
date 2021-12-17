# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:17:48 2019

@author: SONY
"""

# collecting Positive samples: (label as positive in the final model)
import pandas as pd
import numpy as np
pos=pd.read_csv("CosmicMutantExport.tsv",error_bad_lines=False, sep="\t")
pos_arr=np.array(pos)
c=0
for i in range(len(pos_arr)):
     if pos_arr[i,22]==38.0 and ">" in pos_arr[i,17] and pos_arr[i,25] != "y" and pos_arr[i,29]=="Confirmed somatic variant" and "Substitution" in pos_arr[i,19] : #Important!!!
         c=c+1
fin_pos=np.array(np.random.rand (c,35),dtype='str')
j=0
for i in range(len(pos_arr)):
     if  pos_arr[i,22]==38.0 and ">" in pos_arr[i,17] and pos_arr[i,25] !="y" and pos_arr[i,29]=="Confirmed somatic variant" and "Substitution" in pos_arr[i,19]   :
         fin_pos[j,:]=pos_arr[i,:]
         j=j+1
fin_pos_df=pd.DataFrame (fin_pos, columns=pos.columns)
fin_pos_df ['freq_all cancer types'] = fin_pos_df.groupby ('Mutation ID')['Mutation ID'].transform('size')
fin_pos_df=fin_pos_df.drop_duplicates(subset=['Mutation CDS','Mutation genome position','Primary site'] , keep = 'first' , inplace=False)
ID_dict={}
for idd in set(fin_pos_df['Mutation ID']):
                ID_dict[idd]=0
cancertype=list(set(fin_pos_df["Primary site"]))
fin_pos_arr=np.array(fin_pos_df)
for j in cancertype :
      for i in range(len(fin_pos_arr)):
             ID = fin_pos_arr [i,16]
             cancer=fin_pos_arr[i,7]
             if j==cancer:
                 ID_dict[ID]=ID_dict[ID]+1   


#COUNTING FOR GRAPHS#
rec=list(range(1,16))
allcancer_pos_count=[] # for positive samples frequncy in all cancer types
for r in rec:
   c=0
   for i in fin_pos_arr[:,-1]:
         if i>=r:
             c=c+1
   allcancer_pos_count.append(c) 
NumberOfCancers_pos=[] 
for r in rec:
   c=0
   for i in ID_dict.keys():
         if ID_dict[i ]>=r:
             c=c+1
   NumberOfCancers_pos.append(c) 
print(allcancer_pos_count)
print(NumberOfCancers_pos)
# when you decided about which r to choose, then you should select the ones bigger than your
# desired r and put them in a list and then go through the id column and see if the id not in
# your list then drop it or add it to a newly made array!
head=fin_pos_df.columns
w=open("cosmic_coding_all cancer types.csv","w")
w.write(",".join(map(str,head))+'\n')
for i in range(len(fin_pos_arr)):
      if fin_pos_arr[i,-1 ]>=5: #  r=5
            line=fin_pos_arr[i]
            newline=",".join(map(str,line))+'\n'
            w.write(newline)
w.close()
acceptedIDs=[] # list of ids with selected r
for i in ID_dict.keys():
         if ID_dict[i ]>=3: ### r=3
              acceptedIDs.append(i) 
w=open("cosmic_coding_number of cancer types.csv","w")
for j in acceptedIDs:
      for i in range(len(fin_pos_arr)):
               if j == fin_pos_arr [i,16] : 
                    line=fin_pos_arr [i]
                    newline=",".join(map(str,line))+'\n'
                    w.write(newline)
                    print (i)
                    break
w.colse()
