#!/usr/bin/env python
# coding: utf-8

# ***LOADING LIBRARIES***

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


# # IMPORTING LIBRARIES

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


# In[9]:


import os
os.chdir('/home/sanjay/WGA/test/')


# In[10]:


scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
scv.set_figure_params()


# In[12]:


#reading files
file = "scvelo.loom"
adata = scv.read_loom(file)
adata


# In[13]:


#making variable names unique
adata.var_names_make_unique()
scv.pl.proportions(adata)


# # sub-cluster the OSNs

# In[14]:


#reading mature barcodes
m = pd.read_csv('wild_mature.csv')
m = m['Unnamed: 0']
m = m.values.tolist()
a = adata.obs.index
a = a.to_list()


# In[15]:


#intersecting mIOSNs with index of data
l = []
for i in range(len(a)):
    for j in range(len(m)):
        if (a[i] == (m[j])):
            l.append(a[i])


# In[16]:


len(l)


# In[17]:


#reading immature barcodes
im = pd.read_csv('wild_immature.csv')
im = im['Unnamed: 0']
im = im.values.tolist()


# In[18]:


#intersecting iIOSNs with index of data
li = []
for i in range(len(a)):
    for j in range(len(im)):
        if (a[i] == (im[j])):
            li.append(a[i])


# In[19]:


len(li)


# **The selected cell data is retrieved from metadata**

# In[20]:


matcell = l
immatcell = li


# In[21]:


mature = adata[matcell,:]
immature = adata[immatcell, :]


# In[22]:


mature.obs['cluster'] = 'mature'
immature.obs['cluster'] = 'immature'


# In[23]:


final = mature.concatenate(immature)


# In[24]:


i = final.obs.index.tolist()
it = [sub[ : -2] for sub in i]
final.obs.index = it


# # SAVING RAW FILE

# In[25]:


#file with all genes
bdata = final
bdata.write_loom("raw_data_cluster.loom", write_obsm_varm=True)


# In[26]:


final


# # STEP 3: The OSN data being processed for analysis

# In[27]:


#read loom file
f = scv.read_loom("raw_data_cluster.loom")
f


# In[28]:


#make a copy of file with all genes
all_genes = f.copy()
all_genes


# In[29]:


# filter cells and genes
scp.pp.filter_cells(f, min_genes=200)
scp.pp.filter_genes(f, min_cells=3)


# In[30]:


# annotate the group of mitochondrial genes as 'mt'
f.var['mt'] = f.var_names.str.startswith('MT-')  
scp.pp.calculate_qc_metrics(f, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[31]:


#plot total counts, genes_by_counts, and mitochondrial genes
scp.pl.violin(f, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = 'cell_counts.pdf')


# In[32]:


#scatter plot of n genes by count to total count
scp.pl.scatter(f, x='total_counts', y='n_genes_by_counts')


# In[33]:


#filter cells on basis of mitochrial genes and n_genes_by_counts
f = f[f.obs.n_genes_by_counts < 2500, :]
f = f[f.obs.pct_counts_mt < 5, :]
f


# In[118]:


#adding a column of gene and cell ID
f.obs['CellID'] = f.obs.index
f.var['Gene'] = f.var.index


# # SAVING FILTERED FILE

# In[35]:


#writing file with filtered genes
cdata = f
cdata.write_loom("filtered_data_cluster.loom", write_obsm_varm=True)


# In[36]:


#basic preprocessing of data 
#filtering and normalizing data on top 2000 genes, calculate neighbors and moments
scv.pp.filter_and_normalize(f, n_top_genes=2000)
scv.pp.neighbors(f)
scv.pp.moments(f, n_pcs=30, n_neighbors=30)


# In[37]:


#calculate velocity for each cell and velocity graph
scv.tl.velocity(f)
scv.tl.velocity_graph(f)


# In[38]:


# calculating umap coordinates and clustering of cells on the basis of louvain algorithm
scv.tl.umap(f, n_components=2)
scv.tl.louvain(f)


# # READ RAW DATA AND CHECK EXPRESSION

# In[39]:


#reading loom file with filtered genes
ddata = scv.read_loom("filtered_data_cluster.loom", sparse=False)
ddata


# In[40]:


#adding a column to filtered file with louvain cluster number for each cell 
c = f.obs['louvain']
ddata.obs['louvain'] = c


# In[41]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43"]


# In[42]:


#dot plot of immature neuron markers
scp.pl.dotplot(ddata, immature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-immature-exp.pdf')


# In[43]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Myo6", "Stoml3","Gng13","Ech1","Dcx","Gnal"]


# In[44]:


#dotplot of mature neuron markers
scp.pl.dotplot(ddata, mature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-mature-exp.pdf')


# In[45]:


#plotting velocity embedding stream
scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "louvain",save = 'velocity-embedding-stream-steady.pdf')


# In[46]:


#Author provided annotation
scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "cluster",save = 'velocity-embedding-stream-cluster.pdf')


# #Final cluster names
# 

# In[47]:


#naming clusters on the basis of marker genes
f.obs['celltype'] = 'Immature cells'
f.obs['celltype'][f.obs['louvain']=='0'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='4'] = 'Transition cells'
f.obs['celltype'][f.obs['louvain']=='2'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='3'] = 'Mature cells'
f.obs['celltype'][f.obs['louvain']=='1'] = 'Mature cells'


# In[48]:


#Final celltype annotation
scv.pl.velocity_embedding_stream(f, basis = 'umap', color = 'celltype',save= 'velocity_embedding_stream.pdf')


# In[49]:


#adding celltype for filtered gene file
k = f.obs['celltype']
ddata.obs['celltype'] = k


# In[50]:


markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43", "Myo6", "Stoml3",
                           "Gng13","Ech1","Dcx","Gnal"]


# In[51]:


#dot plot of all markers
scp.pl.dotplot(ddata, markers,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot_all_markers.pdf')


# POINT6

# In[52]:


#plotting one exemplary OR
scp.pl.scatter(f, basis = 'umap', color = ['celltype', 'Olfr186'], size = 250, save = 'olfr186-umap.pdf')


# EXTRACTING BARCODES

# In[55]:


trans = f[f.obs['celltype'] == 'Transition cells']
mature = f[f.obs['celltype'] == 'Mature cells']
immature = f[f.obs['celltype'] == 'Immature cells']
trans_bar = pd.DataFrame(trans.obs.index.tolist())
mat_bar = pd.DataFrame(mature.obs.index.tolist())
imm_bar = pd.DataFrame(immature.obs.index.tolist())
trans_bar.to_csv('transition_barcodes.csv')
mat_bar.to_csv('mature_barcodes.csv')
imm_bar.to_csv('immature_barcodeds')


# # CLUSTER SEGREGATION TO FIND THEIR IDENTITY IN R
# 
# 
# 

# IMMATURE

# In[56]:


#subsetting immature cells
clus = (f.obs["celltype"]  == "Immature cells")
clus1 = f[clus].copy()
clus1_ind = clus1.obs_names
clus1_list = clus1_ind.to_list()
clust1 = adata[clus1_list, :]


# Umap coordinates

# In[57]:


#extracting umap coordinates
umap1 = clus1.obsm['X_umap']
umap1 = pd.DataFrame(umap1)
umap1.index = clus1.obs_names
umap1.to_csv('umapimmat.csv')


# TRANSITION

# In[58]:


#subsetting transition cells
clus = (f.obs["celltype"]  == "Transition cells")
clus0 = f[clus].copy()
clus0_ind = clus0.obs_names
clus0_list = clus0_ind.to_list()
#clus0_list = [sub[ : -2] for sub in clus0_ind]
clust0 = adata[clus0_list, :]


# umap

# In[59]:


#extracting umap coordinates
umap0 = clus0.obsm['X_umap']
umap0 = pd.DataFrame(umap0)
umap0.index = clus0.obs_names
umap0.to_csv('umaptrans.csv')


# #cluster mature

# In[60]:


#subsetting mature cells
clus = (f.obs["celltype"]  == "Mature cells")
clus4 = f[clus].copy()
clus4_ind = clus4.obs_names
clus4_list = clus4_ind.to_list()
#clus4_list = [sub[ : -2] for sub in clus4_ind]
clust4 = adata[clus4_list, :]


# umap

# In[61]:


#extracting umap coordinates
umap4 = clus4.obsm['X_umap']
umap4 = pd.DataFrame(umap4)
umap4.index = clus4.obs_names
umap4.to_csv('umapmature.csv')


# # Downstream analysis
# IMPORTANT GENES

# In[62]:


#file with all genes
ddata = f
ddata.write_loom("downstream.loom", write_obsm_varm=True)


# In[63]:


#ranking velocity genes
scv.tl.rank_velocity_genes(f, groupby='celltype', min_corr=.3)
df = scv.DataFrame(f.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv('rank_velocity_genes_wild.csv')


# SPEED AND COHERENCE

# In[64]:


#calculating the rate of differentiation 
scv.tl.velocity_confidence(f)


# In[116]:


keys = 'velocity_length', 'velocity_confidence'
df = f.obs.groupby('celltype')[keys].mean().T


# In[67]:


df.to_csv('velocity_data.csv')


# PAGA VELOCITY GRAPH

# In[117]:


#calculating PAGA for trajectory inference
f.uns['neighbors']['distances'] = f.obsp['distances']
f.uns['neighbors']['connectivities'] = f.obsp['connectivities']
scv.tl.paga(f, groups='celltype')
df = scv.get_df(f, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[69]:


scv.pl.paga(f, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save = 'PAGA.pdf')


# # Dynamical modeling

# In[70]:


#recovering dynamics of the data
scv.tl.recover_dynamics(f)


# In[71]:


#calculating velocity based on dynamical model
scv.tl.velocity(f, mode='dynamical')
scv.tl.velocity_graph(f)


# In[72]:


scv.pl.velocity_embedding_stream(f, basis='umap', color = 'celltype', save = 'differential_velocity_stream.pdf')


# CLUSTER SPECIFIC TOP LIKELIHOOD GENES

# In[73]:


#ranking dynamical genes
scv.tl.rank_dynamical_genes(f, groupby='celltype')
df = scv.get_df(f, 'rank_dynamical_genes/names')
df.head(5)
df.to_csv('rank_dynamical_genes.csv')


# # CELL RANK

# In[74]:


#calculating terminal states
cr.tl.terminal_states(f, weight_connectivities=0.2,cluster_key='celltype')


# In[75]:


cr.tl.initial_states(f, cluster_key='celltype')
#calculate lineages
cr.tl.lineages(f)


# In[76]:


scv.tl.latent_time(f, root_key= 'initial_state_probs', end_key='terminal_states_probs')


# In[77]:


scv.tl.paga(f, groups='celltype', root_key='initial_state_probs', end_key='terminal_states_probs',
            use_time_prior='velocity_pseudotime')


# In[78]:


cr.pl.cluster_fates(f, mode="paga_pie", cluster_key="celltype", basis='umap', 
                    legend_kwargs={'loc': 'top right out'}, legend_loc='top left out', 
                    node_size_scale=5, edge_width_scale=1, max_edge_width=4, title='directed PAGA', save = 'Cellrank-Differential_PAGA.pdf')


# In[79]:


root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]


# In[80]:


f.uns['iroot'] = root_idx


# In[81]:


scp.tl.dpt(f)


# In[82]:


model = cr.ul.models.GAM(f)


# ADVANCE

# In[93]:


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(f)


# In[99]:


#computing transition matrix
vk.compute_transition_matrix()


# In[100]:


from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(f).compute_transition_matrix()


# In[102]:


combined_kernel = 0.8 * vk + 0.2 * ck
root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]
f.uns['iroot'] = root_idx
from cellrank.tl.kernels import PalantirKernel
pk = PalantirKernel(f, time_key='dpt_pseudotime').compute_transition_matrix()


# In[103]:


from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)


# In[105]:


g.compute_schur(n_components=20)
g.compute_macrostates(n_states=2, cluster_key='celltype')
g.set_terminal_states_from_macrostates()
g.compute_absorption_probabilities()
cr.pl.cluster_fates(f, cluster_key='celltype', mode='heatmap', save = 'cellrank-heatmap-clusterfate.pdf')


# # Differential kinetics

# In[106]:


#running differential kinetic test
scv.tl.differential_kinetic_test(f, groupby='celltype')


# TOP LIKELIHOOD genes

# In[107]:


#recovering dynamics
scv.tl.recover_dynamics(f)


# In[108]:


#running differential model
top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(f, groupby='celltype')


# RECOMPUTING VELOCITIES

# In[109]:


#recomputing velocities on basis of differential model
scv.tl.velocity(f, diff_kinetics=True, groupby='celltype')
scv.tl.velocity_graph(f)


# Differential gene expression

# In[110]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=5) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'diffexp5genes.pdf')


# In[111]:


im = ['Tubb5', 'Sox11', 'Tmsb10', 'Ptma', 'Stmn2']
scp.pl.dotplot(f, im,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-immature-diffexp.pdf')


# In[112]:


tr = ['Stmn4', 'Tmsb4x', 'Ptn', 'Gap43', 'Calm2']
scp.pl.dotplot(f, tr,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-trans-diffexp.pdf')


# In[113]:


m = ['S100a5', 'Pcp4l1', 'Gnal', 'Ckb', 'Malat1']
scp.pl.dotplot(f, m,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-mature-diffexp.pdf')


# velocity length and confidence

# In[114]:


vl = f.obs['velocity_length']
vl.to_csv('velocity_length_wild.csv')


# In[115]:


vc = f.obs['velocity_confidence']
vc.to_csv('velocity_confi_wild.csv')


# In[ ]:




