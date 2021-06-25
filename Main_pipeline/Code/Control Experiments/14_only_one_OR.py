#!/usr/bin/env python
# coding: utf-8

# In[49]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


# In[50]:


os.chdir('/home/sidrah19220/nmd_review/reviewer2/one_OR')


# # RUN THE EDI 

# In[51]:


eud_orig = pd.read_csv('matrix-1-euclidean-correlation.csv', index_col =0)
exp = pd.read_csv('matrix-5-Expression.csv', index_col = 0)


# In[52]:


eud = eud_orig
eu = eud_orig


# In[53]:


ors=list(exp.columns.values.tolist())
cells=exp.index.values
#exp.index = cells
tr_exp=exp.T


# In[54]:


all_ors={}
zero_exp=[]
for cell in cells:
  all_ors[cell]=[]
  cell_exps=sorted(tr_exp[cell].tolist(), reverse=True)
  if cell_exps[0]>0:
    if tr_exp[cell].tolist().count(max(cell_exps)) > 8:
      for i in range(tr_exp[cell].tolist().count(max(tr_exp[cell]))):
        all_ors[cell].append(tr_exp.sort_values(by=cell, ascending=False).index[i])
    else:
      for i in range(8):
        if cell_exps[i]>0:
          all_ors[cell].append(tr_exp.sort_values(by=cell, ascending=False).index[i])
  else:
    zero_exp.append(cell)


# In[55]:


cell_imt={}
cell_imm={}
cell_mat={}
for key,value in all_ors.items():
  if key.startswith("IMT"):
    cell_imt[key]=[]
    for val in value:
      if val.startswith("Olfm"):
        continue
      else:
        cell_imt[key].append(val)
  if key.startswith("MA"):
    cell_mat[key]=[]
    for val in value:
      if val.startswith("Olfm"):
        continue
      else:
        cell_mat[key].append(val)
  if key.startswith("IMM"):
    for val in value:
      cell_imm[key]=[]
      if val.startswith("Olfm"):
        continue
      else:
        cell_imm[key].append(val)


# In[56]:


cell_imt['IMT CELL280']


# In[57]:


cell_imm_imt={}
for key,value in all_ors.items():
  if key.startswith("IM"):
    cell_imm_imt[key]=[]
    for val in value:
      if val.startswith("Olfm"):
        continue
      else:
        cell_imm_imt[key].append(val)


# In[58]:


di3={}
for key,value in cell_imm_imt.items():
    di4={}
    if len(value)<1:
        continue
    elif len(value)>1:
        continue
    else:
        for val in value:
            di4[val]=tr_exp[key][val]
            di3[key]=di4


# In[59]:


cell_imm_imt = di3


# In[60]:


di1={}
for key,value in cell_imm_imt.items():
  di2={}
  if len(value)<1:
    continue
  else:
    for val in value:
      di2[val]=tr_exp[key][val]
    di1[key]=di2


# In[61]:


cell_imm_imt


# In[62]:


#GENERATE UNIQUE HITS
os.chdir('/home/sidrah19220/nmd_review/reviewer2/one_OR')
ratio_data={}
unique_number = {}
unique_percentage = {}
for folder in range(1):
    di='test_'+str(folder)
    os.mkdir(di)
    os.chdir(di) 
        #A FILE IS GENERATED WITH ALL INFORMATION OF mOSN, iOSN/tOSN, THE EUCLIDEAN DISTANCE AND OR EXPRESSED
    with open("Out_data_eud_"+str(folder)+".csv","w")as fout:
        fout.write("Olfactory Receptor,Mature Cell ID,Immature/Transition Cell ID, ED,Dist. Rank,Receptor Expression Rank,Receptor Expression,Expression Diffrence from rank1\n")
        empty_mat=[]
        mat_cell=eud.index
        imm_imt_cells=list(eud.columns.values.tolist())
        for i in range(len(mat_cell)):
            mat1={}
            for j in range(len(imm_imt_cells)):
                mat1[eud.iloc[i][j]]=imm_imt_cells[j]
            euds=sorted(mat1.keys(), reverse=False)
            for j in range(len(euds)):
                if len(cell_mat[mat_cell[i]])==0:
                    empty_mat.append(mat_cell[i])
                else:
                    if mat1[euds[j]] in cell_imm_imt.keys():
                        if cell_mat[mat_cell[i]][0] in cell_imm_imt[mat1[euds[j]]]:
                            for k in range(len(cell_imm_imt[mat1[euds[j]]])):
                                if list(cell_imm_imt[mat1[euds[j]]].keys())[k]==cell_mat[mat_cell[i]][0]:
                                    fout.write(str(list(cell_imm_imt[mat1[euds[j]]].keys())[k])+','+str(mat_cell[i])+','+str(mat1[euds[j]])+','+str(euds[j])+','+str(j)+','+str(k)+','+str(di1[mat1[euds[j]]][list(cell_imm_imt[mat1[euds[j]]].keys())[k]])+','+str((di1[mat1[euds[j]]][list(cell_imm_imt[mat1[euds[j]]].keys())[0]])-di1[mat1[euds[j]]][list(cell_imm_imt[mat1[euds[j]]].keys())[k]])+'\n')

        
        #CALCULATING PERCENTILE AT DIFFERENT THRESHOLD FOR EUC DISTANCE
    eu = eud
    li = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for ki in li:
        q = eu.quantile(ki , axis = 1)
        q1 = q.tolist()
        qq = pd.DataFrame(q)

        idd = eu.index.tolist()
        ccl = eu.columns.tolist()
        d4=dict.fromkeys(idd)

        #loop to generate imm cells
        k = []
        for i in range(len(idd)):
            for j in range(len(ccl)):
                if eu.iloc[i][j] < q1[i]:
                    k.append(ccl[j])
            d4[idd[i]] = k
            k=[]

        p = []
        for i in idd:
            p.append(len(d4[i]))
        cells = pd.DataFrame.from_dict(d4, orient='index')
        cells.to_csv(str(ki)+'percentile_euc_total_immature.csv')
    
    
    #ASSIGNING THE STATUS POSITIVE AND NEGATIVE ON BASIS OF THE mOSN THAT LIE IN THAT PARTICULAR THRESHOLD
    for ki in li:
        p=pd.read_csv(str(ki)+'percentile_euc_total_immature.csv', index_col = 0)
        o = pd.read_csv('Out_data_eud_'+str(folder)+'.csv')
        o = o[o[' ED']!=0]
        dictt={}
        v =  o['Mature Cell ID'].tolist()
        v2 = o['Immature/Transition Cell ID'].tolist()
        for indx in range(len(v)):
            if v[indx]==v[indx-1]:
                dictt[v[indx]].append(v2[indx])
            else:
                dictt[v[indx]]=[]
                dictt[v[indx]].append(v2[indx])
        ke=[]
        for key in dictt:
            ke.append(key)
        pp = p.loc[ke].T
        ppp1 = pp.to_dict('list')
        dup = o
        ind = dup.index.tolist()
        stat=[]
        for i in range(0,len(dup['Immature/Transition Cell ID'])):
            if dup['Immature/Transition Cell ID'][ind[i]] in ppp1[dup['Mature Cell ID'][ind[i]]]:
                stat.append('Positive')
            else:
                stat.append('Negative')

        dup['status'] = stat
        #dup.to_csv('Final_spl_10q.csv')
        dup.to_csv(str(ki)+'percentile_neg_pos_out_data.csv')
        euk = dup[(dup['status']=='Positive') & (dup['Expression Diffrence from rank1']==0)]
        euk.to_csv(str(ki)+'_final_res.csv')

        #CREATING A DATAFRAME AT DIFFERENT PERCENTILE SHOWING RATIO OF EXPRESSED ORs HAVING HIGHEST EXPRESSION TO NON-EXPRESSED
    DF4=pd.DataFrame()
    for ki in li:
        x = pd.read_csv(str(ki)+'_final_res.csv', index_col = 0)
        o = pd.read_csv('Out_data_eud_'+str(folder)+".csv", index_col =0)
        o = o[o[' ED']!=0]
        eu = eud
        mat_uniq = eu.index.tolist()
        dictt={}
        v =  x['Mature Cell ID'].tolist()
        v2 = x['Immature/Transition Cell ID'].tolist()

        for indx in range(len(v)):
            if v[indx]==v[indx-1]:
                if v[indx] not in dictt.keys():
                    dictt[v[indx]]=[]
                dictt[v[indx]].append(v2[indx])
            else:
                dictt[v[indx]]=[]
                dictt[v[indx]].append(v2[indx])
        xa_dict={}
        for ma in mat_uniq:
            if ma in dictt.keys():
                xa_dict[ma]=len(dictt[ma])
            else:
                xa_dict[ma]=0
        df2 = pd.DataFrame.from_dict(xa_dict, orient='index')
        print(df2.value_counts())

        pos_neg = pd.read_csv(str(ki)+'percentile_neg_pos_out_data.csv')
        pos_neg2 = pos_neg[pos_neg['status'] == 'Positive']
        dictt2={}
        v =  pos_neg2['Mature Cell ID'].tolist()
        v2 = pos_neg2['Immature/Transition Cell ID'].tolist()

        for indx in range(len(v)):
            if v[indx]==v[indx-1]:
                if v[indx] not in dictt2.keys():
                    dictt2[v[indx]]=[]
                dictt2[v[indx]].append(v2[indx])
            else:
                dictt2[v[indx]]=[]
                dictt2[v[indx]].append(v2[indx])

        x_dict={}
        for ma in mat_uniq:
            if ma in dictt2.keys():
                x_dict[ma]=len(dictt2[ma])
            else:
                x_dict[ma]=0

        df3 = pd.DataFrame.from_dict(x_dict, orient='index')
        print(df3.value_counts())
        df3.columns = ['x']
        percentile = pd.read_csv(str(ki)+'percentile_euc_total_immature.csv', index_col =0)
        xny=len(percentile.columns.values.tolist())
        print(xny)
        if len(percentile.index.values.tolist())==len(df3):
            df3['y'] = xny - df3['x'] 
        xa = df2[0].values.tolist()
        df3['xa'] = xa
        df3['xay'] = df3['xa']/df3['y']
        df3.to_csv(str(ki*100)+'percentile_boxplot.csv')

        DF4[str(ki*100)+'p']=df3['xay']

    #DF4.to_csv('plot_boxplot_euclidean.csv')
    p  = o.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')['Mature Cell ID'].tolist()
    DF4.loc[p].boxplot()
    plt.xlabel('percentile')
    plt.ylabel('Xa/y')
    plt.title('Normalised data')
    plt.savefig('normalised_data_euclidean.pdf')
    plt.show()

    o_unique = o.drop_duplicates(subset = 'Mature Cell ID', keep ='first')
    mat_uniq= o_unique['Mature Cell ID'].values.tolist()
    dictt={}
    v =  o['Mature Cell ID'].tolist()
    v2 = o['Immature/Transition Cell ID'].tolist()

    for indx in range(len(v)):
        if v[indx]==v[indx-1]:

            dictt[v[indx]].append(v2[indx])
        else:
            dictt[v[indx]]=[]
            dictt[v[indx]].append(v2[indx])

    f1 = pd.read_csv('0.01_final_res.csv')
    f2 = pd.read_csv('0.02_final_res.csv')
    f5 = pd.read_csv('0.05_final_res.csv')
    f10 = pd.read_csv('0.1_final_res.csv')
    f20 = pd.read_csv('0.2_final_res.csv')
    f30 = pd.read_csv('0.3_final_res.csv')
    f50 = pd.read_csv('0.5_final_res.csv')
    f60 = pd.read_csv('0.6_final_res.csv')
    f70 = pd.read_csv('0.7_final_res.csv')
    f80 = pd.read_csv('0.8_final_res.csv')
    f90 = pd.read_csv('0.9_final_res.csv')
    f100 = pd.read_csv('1_final_res.csv')

        #PLOT OF NUMBER OF HITS
    
    liii=[f1, f2, f5, f10, f20, f30, f50, f60, f70, f80, f90, f100]
    valss=[]
    for i in liii:
        valss.append(len(i))

    x = ['1p', '2p', '5p', '10p', '20p', '30p', '50p', '60p', '70p', '80p', '90p', '100p']
    plt.bar(x,valss)
    plt.title('Count of hits')
    plt.xlabel('percentile')
    plt.ylabel('# of cells')
    plt.savefig('number of hits euclidean.pdf')
    plt.show()

    unique_f1 = f1.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f2 = f2.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f5 = f5.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f10 = f10.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f20 = f20.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f30 = f30.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f50 = f50.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f60 = f60.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f70 = f70.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f80 = f80.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f90 = f90.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')
    unique_f100 = f100.drop_duplicates(subset = 'Mature Cell ID', keep = 'first')

    
    #PLOT OF NUMBER OF UNIQUE mOSNs 

    liii=[unique_f1, unique_f2, unique_f5, unique_f10, unique_f20, unique_f30, unique_f50, unique_f60, unique_f70, unique_f80,
          unique_f90, unique_f100]
    valss_unique=[]
    for i in liii:
        valss_unique.append(len(i))

    unique_number[str(folder)]=valss_unique
    
    gra = pd.DataFrame(valss_unique)
    gra.to_csv('unique-count.csv')
    
    
    x = ['1p', '2p', '5p', '10p', '20p', '30p', '50p', '60p', '70p', '80p', '90p', '100p']
    plt.bar(x,valss_unique)
    plt.title('Count of unique hits')
    plt.xlabel('percentile')
    plt.ylabel('# of cells')
    plt.savefig('number of unique hits euclidean.pdf')
    plt.show()
    
    
    #PERCENTAGE OF UNIQUE mOSNs
    valssf=[]
    for v in valss_unique:
        valssf.append((v/1559)*100)
        
    unique_percentage[str(folder)]=valssf
    
    
    gra = pd.DataFrame(valssf)
    gra.to_csv('unique-percentage.csv')
    

    plt.bar(x,valssf)
    plt.title('Percentage of unique counts')
    plt.xlabel('Percentile')
    plt.ylabel('Percentage')
    plt.savefig('percentage')
    plt.yticks([0, 20, 40, 60, 80, 100])
    plt.savefig('percentage_test_divbyunique.pdf')
    plt.show()
    
    os.chdir('/home/sidrah19220/nmd_review/reviewer2/one_OR')


# In[46]:


list(cell_imm_imt[mat1[euds[j]]].keys())[0],cell_mat[mat_cell[i]][0]


# In[43]:


mat1


# In[ ]:





# In[19]:


one_bar = []
two_bar =[]

for i in range(len(one)):
    one_bar.append(bar.loc[one[i]])
for j in range(len(two)):
    two_bar.append(bar.loc[two[j]])


# In[20]:


one_barr = pd.DataFrame(one_bar)
two_barr = pd.DataFrame(two_bar)
one_barr.to_csv('One_barr.csv')
two_barr.to_csv('two_barr.csv')


# READ LOOM FILE AND RUN SCVELO TO EXTRACT VELOCITY VECTORS

# In[22]:


import scvelo as scv
import pandas as pd
#import velocyto as vcy
import time
import scipy
import pickle
import numpy as np
import pandas as pd
import scanpy as scp
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib_venn as ven
from matplotlib.gridspec import GridSpec

import sklearn
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error, accuracy_score, cohen_kappa_score, confusion_matrix, recall_score, classification_report, precision_score
from scipy.stats import spearmanr, zscore
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import Lasso, Ridge, ElasticNet
import scvelo as scv
scv.logging.print_version()
import jgraph as ig
from sklearn.svm import SVC


# In[23]:


#reading files
file = "wild_type.loom"
adata = scv.read_loom(file)
adata


# In[24]:


#making variable names unique
adata.var_names_make_unique()
scv.pl.proportions(adata)


# # sub-cluster the OSNs

# In[25]:


#reading mature barcodes
m = pd.read_csv('wild_mature.csv')
m = m['Unnamed: 0']
m = m.values.tolist()
a = adata.obs.index
a = a.to_list()


# In[26]:


#intersecting mIOSNs with index of data
l = []
for i in range(len(a)):
    for j in range(len(m)):
        if (a[i] == (m[j])):
            l.append(a[i])


# In[27]:


#reading immature barcodes
im = pd.read_csv('wild_immature.csv')
im = im['Unnamed: 0']
im = im.values.tolist()


# In[28]:


#intersecting iIOSNs with index of data
li = []
for i in range(len(a)):
    for j in range(len(im)):
        if (a[i] == (im[j])):
            li.append(a[i])


# In[29]:


matcell = l
immatcell = li


# In[30]:


mature = adata[matcell,:]
immature = adata[immatcell, :]


# In[31]:


mature.obs['cluster'] = 'mature'
immature.obs['cluster'] = 'immature'


# In[32]:


final = mature.concatenate(immature)
i = final.obs.index.tolist()
it = [sub[ : -2] for sub in i]
final.obs.index = it


# In[33]:


#read loom file
f = scv.read_loom("raw_data_cluster.loom")
f


# # DATA PREPROCESSING

# In[35]:


f = final


# In[36]:


#make a copy of file with all genes
all_genes = f.copy()
all_genes


# In[37]:


# filter cells and genes
scp.pp.filter_cells(f, min_genes=200)
scp.pp.filter_genes(f, min_cells=3)


# In[38]:


# annotate the group of mitochondrial genes as 'mt'
f.var['mt'] = f.var_names.str.startswith('MT-')  
scp.pp.calculate_qc_metrics(f, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[39]:


#plot total counts, genes_by_counts, and mitochondrial genes
scp.pl.violin(f, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = 'cell_counts.pdf')


# In[40]:


#scatter plot of n genes by count to total count
scp.pl.scatter(f, x='total_counts', y='n_genes_by_counts')


# In[41]:


#filter cells on basis of mitochrial genes and n_genes_by_counts
f = f[f.obs.n_genes_by_counts < 2500, :]
f = f[f.obs.pct_counts_mt < 5, :]
f


# In[42]:


#adding a column of gene and cell ID
print(f.obs.index)
f.obs['CellID'] = f.obs.index
print(f.obs)
print(f.var.index)
f.var['Gene'] = f.var.index
print(f.var)


# In[43]:


#writing file with filtered genes
cdata = f
cdata.write_loom("filtered_data_cluster.loom", write_obsm_varm=True)


# # SCVELO

# In[44]:


#basic preprocessing of data 
#filtering and normalizing data on top 2000 genes, calculate neighbors and moments
scv.pp.filter_and_normalize(f, n_top_genes=2000)
scv.pp.neighbors(f)
scv.pp.moments(f, n_pcs=30, n_neighbors=30)


# In[45]:


#calculate velocity for each cell and velocity graph
scv.tl.velocity(f)
scv.tl.velocity_graph(f)


# In[46]:


# calculating umap coordinates and clustering of cells on the basis of louvain algorithm
scv.tl.umap(f, n_components=2)
scv.tl.louvain(f)


# In[47]:


#reading loom file with filtered genes
ddata = scv.read_loom("filtered_data_cluster.loom", sparse=False)
ddata


# In[48]:


#adding a column to filtered file with louvain cluster number for each cell 
c = f.obs['louvain']
ddata.obs['louvain'] = c


# In[49]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43"]
scp.pl.violin(ddata, keys=immature_neuron_markers, groupby='louvain', save = 'immature_neurons_exp.pdf')


# In[50]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Myo6", "Stoml3","Gng13","Ech1","Dcx","Gnal"]
scp.pl.violin(ddata, keys=mature_neuron_markers, groupby='louvain', save = 'mature-exp.pdf')


# In[51]:


#velocity embedding plotting
scv.pl.velocity_embedding(f, basis = 'umap', color = 'louvain', save = 'velocity-embedding.pdf')


# In[52]:


#plotting velocity embedding stream
scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "louvain",save = 'velocity-embedding-stream-steady.pdf')


# In[53]:


#naming clusters on the basis of marker genes
f.obs['celltype'] = 'Immature cells'
f.obs['celltype'][f.obs['louvain']=='0'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='4'] = 'Transition cells'
f.obs['celltype'][f.obs['louvain']=='2'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='3'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='1'] = 'Mature cells'


# In[54]:


scv.pl.velocity_embedding_stream(f, basis = 'umap', color = 'celltype',save= 'velocity_embedding_stream.pdf')


# In[55]:


#adding celltype for filtered gene file
k = f.obs['celltype']
ddata.obs['celltype'] = k


# In[57]:


bar_one = one_barr['0'].values.tolist()
bar_two = two_barr['0'].values.tolist()


# In[59]:


imm_one = f[bar_one, :]
imm_two = f[bar_two, :]


# In[72]:


l1 = imm_one.obs_names.tolist()
l2 = imm_two.obs_names.tolist()


# In[73]:


imm_one_bar = np.array(imm_one.layers['velocity'])
imm_two_bar = np.array(imm_two.layers['velocity'])


# In[74]:


vel1 = pd.DataFrame(imm_one_bar)
vel2 = pd.DataFrame(imm_two_bar)


# In[86]:


g1 = imm_one.var_names.tolist()
g2 = imm_two.var_names.tolist()


# In[78]:


vel1.index = l1
vel2.index = l2


# In[87]:


vel1.columns = g1
vel2.columns = g2


# In[89]:


vel1.to_csv('immature_one_or_vel_vector.csv')
vel2.to_csv('immature_two_or_vel_vector.csv')


# In[136]:


#OR list saved

pd.DataFrame(exp.columns.tolist()).to_csv('OR_list.csv')

