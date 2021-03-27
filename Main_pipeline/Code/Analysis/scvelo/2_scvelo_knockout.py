#!/usr/bin/env python
# coding: utf-8

# # ***LOADING LIBRARIES***

# In[1]:


pip install cellrank --quiet


# In[2]:


get_ipython().system('pip install scvelo --upgrade --quiet')


# In[3]:


pip install python-igraph


# In[4]:


pip install jgraph --upgrade --quiet


# In[5]:


pip install louvain


# #IMPORTING LIBRARIES

# In[6]:


import cellrank as cr


# In[7]:


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


# In[8]:


import os
os.chdir('/home/sanjay/WGA/test/test_knock')


# In[9]:


scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
scv.set_figure_params()


# In[10]:


#reading loom file
file = "knock_type.loom"
adata = scv.read_loom(file)
adata


# In[11]:


#making variable names unique
adata.var_names_make_unique()
scv.pl.proportions(adata)


# # sub-cluster the OSNs

# In[12]:


#reading mature barcodes
m = pd.read_csv('knock_mature.csv')
m = m['Unnamed: 0']
m = m.values.tolist()
a = adata.obs.index
a = a.to_list()


# In[13]:


#intersecting mature barcodes with data 
l = []
for i in range(len(a)):
    for j in range(len(m)):
        if (a[i] == (m[j])):
            l.append(a[i])


# In[14]:


len(l)


# In[15]:


#reading immature barcodes
im = pd.read_csv('knock_immature.csv')
im = im['Unnamed: 0']
im = im.values.tolist()


# In[16]:


#intersecting immature barcodes with data
li = []
for i in range(len(a)):
    for j in range(len(im)):
        if (a[i] == (im[j])):
            li.append(a[i])


# In[17]:


len(li)


# **The selected cell data is retrieved from metadata**

# In[18]:


matcell = l
immatcell = li


# In[19]:


mature = adata[matcell,:]
immature = adata[immatcell, :]


# In[20]:


mature.obs['cluster'] = 'mature'
immature.obs['cluster'] = 'immature'


# In[21]:


final = mature.concatenate(immature)


# In[22]:


i = final.obs.index.tolist()
it = [sub[ : -2] for sub in i]
final.obs.index = it


# # SAVING RAW FILE

# In[23]:


#file with all genes
bdata = final
bdata.write_loom("raw_data_cluster.loom", write_obsm_varm=True)


# In[24]:


final


# # STEP 3: The OSN data being processed for analysis

# In[25]:


#reading file with all genes
f = scv.read_loom("raw_data_cluster.loom")
f


# In[26]:


# filtering cells and genes
scp.pp.filter_cells(f, min_genes=200)
scp.pp.filter_genes(f, min_cells=3)


# In[27]:


# annotate the group of mitochondrial genes as 'mt'
f.var['mt'] = f.var_names.str.startswith('MT-')  
scp.pp.calculate_qc_metrics(f, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[28]:


#filtering cells on basis of mitochondrial genes and n genes by counts
f = f[f.obs.n_genes_by_counts < 2500, :]
f = f[f.obs.pct_counts_mt < 5, :]
f


# In[29]:


#adding column of cellID and gene
f.obs['CellID'] = f.obs.index
f.var['Gene'] = f.var.index


# # SAVING FILTERED FILE

# In[30]:


#file with filtered genes
cdata = f
cdata.write_loom("filtered_data_cluster.loom", write_obsm_varm=True)


# In[31]:


#basic preprocessing of data
#filtering and normalizing data, calculating neighbors and moments
scv.pp.filter_and_normalize(f, n_top_genes=2000)
scv.pp.neighbors(f)
scv.pp.moments(f, n_pcs=30, n_neighbors=30)


# In[32]:


#calculating velocity and velocity graph on basis of sparse matrix with cosine correlation
scv.tl.velocity(f)
scv.tl.velocity_graph(f)


# In[33]:


#calculating umap coordinates and clustering using louvain algorithm
scv.tl.umap(f, n_components=2)
scv.tl.louvain(f)


# # READ RAW DATA AND CHECK EXPRESSION

# In[34]:


#reading loom file with filtered genes
ddata = scv.read_loom("filtered_data_cluster.loom", sparse=False)
ddata


# In[35]:


#adding louvain cluster info to filtered file
c = f.obs['louvain']
ddata.obs['louvain'] = c


# In[36]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43"]


# In[37]:


scp.pl.dotplot(ddata, immature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-immature-exp.pdf')


# In[38]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Myo6", "Stoml3","Gng13","Ech1","Dcx","Gnal"]


# In[39]:


scp.pl.dotplot(ddata, mature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-mature-exp.pdf')


# In[40]:


scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "louvain",save = 'velocity-embedding-stream-steady.pdf')


# In[41]:


scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "cluster",save = 'velocity-embedding-stream-cluster.pdf')


# #Final cluster names
# 

# In[42]:


#naming clusters on basis of expression of markers
f.obs['celltype'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='1'] = 'Immature cells'


# In[43]:


scv.pl.velocity_embedding_stream(f, basis = 'umap', color = 'celltype',save= 'velocity_embedding_stream.pdf')


# In[44]:


#adding celltype info to filtered file data
k = f.obs['celltype']
ddata.obs['celltype'] = k


# In[45]:


markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43", "Myo6", "Stoml3",
                           "Gng13","Ech1","Dcx","Gnal"]


# In[46]:


scp.pl.dotplot(ddata, markers,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot_all_markers.pdf')


# # CLUSTER SEGREGATION TO FIND THEIR IDENTITY IN R
# 
# 
# 

# IMMATURE

# In[47]:


#extracting iOSN cell data
clus = (f.obs["celltype"]  == "Immature cells")
clus1 = f[clus].copy()
clus1_ind = clus1.obs_names
clus1_list = clus1_ind.to_list()
#clus1_list = [sub[ : -2] for sub in clus1_ind]
clust1 = adata[clus1_list, :]


# Umap coordinates

# In[48]:


#extract umap coordinates
umap1 = clus1.obsm['X_umap']
umap1 = pd.DataFrame(umap1)
umap1.index = clus1.obs_names
umap1.to_csv('umapimmat.csv')


# spliced of mature --- 1

# In[49]:


#spliced data
spl1 = clust1.layers['spliced']
spl1 = spl1.todense()
spl1 = pd.DataFrame(spl1)
spl1.columns = clust1.var_names
spl1.index = clust1.obs_names
spl1.to_csv('splimmat.csv')


# Unspliced of mature

# In[50]:


#unspliced data
uspl1 = clust1.layers['unspliced']
uspl1 = uspl1.todense()
uspl1 = pd.DataFrame(uspl1)
uspl1.columns = clust1.var_names
uspl1.index = clust1.obs_names
uspl1.to_csv('usplimma.csv')


# #cluster 4

# In[51]:


#mOSN data
clus = (f.obs["celltype"]  == "Mature cells")
clus4 = f[clus].copy()
clus4_ind = clus4.obs_names
clus4_list = clus4_ind.to_list()
#clus4_list = [sub[ : -2] for sub in clus4_ind]
clust4 = adata[clus4_list, :]


# umap

# In[52]:


#umap coordinates
umap4 = clus4.obsm['X_umap']
umap4 = pd.DataFrame(umap4)
umap4.index = clus4.obs_names
umap4.to_csv('umapmature.csv')


# spliced

# In[53]:


#spliced data
spl4 = clust4.layers['spliced']
spl4 = spl4.todense()
spl4 = pd.DataFrame(spl4)


# In[54]:


spl4.columns = clust4.var_names
spl4.index = clust4.obs_names


# In[55]:


spl4.to_csv('splmature.csv')


# unspliced

# In[56]:


#unspliced data
uspl4 = clust4.layers['unspliced']
uspl4 = uspl4.todense()
uspl4 = pd.DataFrame(uspl4)
uspl4.columns = clust4.var_names
uspl4.index = clust4.obs_names
uspl4.to_csv('usplmature.csv')


# # Downstream analysis
# IMPORTANT GENES

# In[57]:


#calculating rank velocity genes
scv.tl.rank_velocity_genes(f, groupby='celltype', min_corr=.3)
df = scv.DataFrame(f.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv('rank-velocity-genes-knock.csv')


# SPEED AND COHERENCE

# In[58]:


#calculating rate of differentiation
scv.tl.velocity_confidence(f)


# In[59]:


keys = 'velocity_length', 'velocity_confidence'
df = f.obs.groupby('celltype')[keys].mean().T


# In[60]:


df.to_csv('velocity-length-and-conf.csv')


# In[61]:


vl = f.obs['velocity_length']
vl.to_csv('velocity-length-knock.csv')


# In[62]:


vm = f.obs['velocity_confidence']
vm.to_csv('velocity-confi-knock.csv')


# PAGA VELOCITY GRAPH

# In[63]:


#calculate PAGA for trajectory inference
f.uns['neighbors']['distances'] = f.obsp['distances']
f.uns['neighbors']['connectivities'] = f.obsp['connectivities']
scv.tl.paga(f, groups='celltype')
df = scv.get_df(f, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[64]:


scv.pl.paga(f, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save = 'PAGA.pdf')


# # Dynamical modeling

# In[65]:


#recovering dynamics of data
scv.tl.recover_dynamics(f)


# In[66]:


#calculating velocity on dynamical model
scv.tl.velocity(f, mode='dynamical')
scv.tl.velocity_graph(f)


# In[67]:


scv.pl.velocity_embedding_stream(f, basis='umap', color = 'celltype', save = 'differential_velocity_stream.pdf')


# In[68]:


#calculating latent time
scv.tl.latent_time(f)


# CLUSTER SPECIFIC TOP LIKELIHOOD GENES

# In[69]:


#ranking dynamical genes
scv.tl.rank_dynamical_genes(f, groupby='celltype')
df = scv.get_df(f, 'rank_dynamical_genes/names')


# In[70]:


df.to_csv('rank-dynamical-genes-kno.csv')


# # CELL RANK

# In[71]:


#calculate terminal states
cr.tl.terminal_states(f, weight_connectivities=0.2,cluster_key='celltype')


# In[72]:


#calculate initial states
cr.tl.initial_states(f, cluster_key='celltype')
#calculate lineages
cr.tl.lineages(f)


# In[73]:


scv.tl.latent_time(f, root_key= 'initial_state_probs', end_key='terminal_states_probs')


# In[74]:


#calculate paga for trajectory inference using cellrank
scv.tl.paga(f, groups='celltype', root_key='initial_state_probs', end_key='terminal_states_probs',
            use_time_prior='velocity_pseudotime')


# In[75]:


cr.pl.cluster_fates(f, mode="paga_pie", cluster_key="celltype", basis='umap', 
                    legend_kwargs={'loc': 'top right out'}, legend_loc='top left out', 
                    node_size_scale=5, edge_width_scale=1, max_edge_width=4, title='directed PAGA', save = 'Cellrank-Differential_PAGA.pdf')


# In[76]:


root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]


# In[77]:


f.uns['iroot'] = root_idx


# In[78]:


scp.tl.dpt(f)


# In[79]:


model = cr.ul.models.GAM(f)


# ADVANCE

# In[80]:


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(f)


# In[81]:


vk.compute_transition_matrix()


# In[82]:


from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(f).compute_transition_matrix()


# In[83]:


combined_kernel = 0.8 * vk + 0.2 * ck
root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]
f.uns['iroot'] = root_idx


# In[84]:


from cellrank.tl.kernels import PalantirKernel
pk = PalantirKernel(f, time_key='dpt_pseudotime').compute_transition_matrix()


# In[85]:


from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)


# In[86]:


g.compute_schur(n_components=20)


# In[87]:


g.compute_macrostates(n_states=2, cluster_key='celltype')


# In[88]:


g.set_terminal_states_from_macrostates()


# In[89]:


g.compute_absorption_probabilities()


# In[90]:


cr.pl.cluster_fates(f, cluster_key='celltype', mode='heatmap', save = 'cellrank-heatmap-clusterfate.pdf')


# # Differential kinetics

# In[91]:


#running differential kinetic test
scv.tl.differential_kinetic_test(f, groupby='celltype')


# TOP LIKELIHOOD genes

# In[92]:


#recovering dynamics of the data
scv.tl.recover_dynamics(f)


# In[93]:


top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(f, groupby='celltype')


# RECOMPUTING VELOCITIES

# In[94]:


#recomputing velocity using differential model
scv.tl.velocity(f, diff_kinetics=True, groupby='celltype')
scv.tl.velocity_graph(f)


# DIFFERENTIAL GENE EXPRESSION

# In[95]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=5) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'diffexp5genes.pdf')


# In[96]:


m = ['S100a5', 'Stoml3', 'Mt1', 'Cartpt', 'Tmbim6']
scp.pl.dotplot(f, m,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-mature-diffexp.pdf')


# In[97]:


im = ['Tmsb10', 'Tubb5', 'Sox11', 'Ptma', 'Stmn2']
scp.pl.dotplot(f, im,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-immature-diffexp.pdf')


# In[ ]:




