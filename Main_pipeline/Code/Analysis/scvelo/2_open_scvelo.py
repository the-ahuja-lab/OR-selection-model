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

# In[1]:


import cellrank as cr


# In[2]:


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


# In[3]:


import os
os.chdir('/home/sidrah19220/nmd/mouse/analysis/scvelo/open')


# In[4]:


scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
scv.set_figure_params()


# In[5]:


file = "open.loom"
adata = scv.read_loom(file)
adata


# In[6]:


adata.var_names_make_unique()
scv.pl.proportions(adata)


# In[7]:


f = adata


# In[8]:


# filter cells and genes
scp.pp.filter_cells(f, min_genes=200)
scp.pp.filter_genes(f, min_cells=3)


# In[9]:


f.var['mt'] = f.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
scp.pp.calculate_qc_metrics(f, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[10]:


scp.pl.violin(f, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = 'cell_counts.pdf')


# In[11]:


scp.pl.scatter(f, x='total_counts', y='n_genes_by_counts')


# In[12]:


f = f[f.obs.n_genes_by_counts < 2500, :]
f = f[f.obs.pct_counts_mt < 5, :]
f


# In[13]:


print(f.obs.index)
f.obs['CellID'] = f.obs.index
print(f.obs)
print(f.var.index)
f.var['Gene'] = f.var.index
print(f.var)


# # SAVING FILTERED FILE

# In[14]:


#file with all genes
cdata = f
cdata.write_loom("filtered_data_cluster.loom", write_obsm_varm=True)


# In[15]:


#basic preprocessing of data
scv.pp.filter_and_normalize(f, n_top_genes=2000)
scv.pp.neighbors(f)
scv.pp.moments(f, n_pcs=30, n_neighbors=30)


# In[16]:


scv.tl.velocity(f)
scv.tl.velocity_graph(f)


# In[17]:


scv.tl.umap(f, n_components=2)
scv.tl.louvain(f)


# # READ RAW DATA AND CHECK EXPRESSION

# In[18]:


ddata = scv.read_loom("filtered_data_cluster.loom", sparse=False)
ddata


# In[19]:


c = f.obs['louvain']
ddata.obs['louvain'] = c


# In[20]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Gap43"]
scp.pl.violin(ddata, keys=immature_neuron_markers, groupby='louvain', save = 'immature_neurons_exp.pdf')


# In[21]:


scp.pl.dotplot(ddata, immature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-immature-exp.pdf')


# In[22]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Omp"]
scp.pl.violin(ddata, keys=mature_neuron_markers, groupby='louvain', save = 'mature-exp.pdf')


# In[23]:


scp.pl.dotplot(ddata, mature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvain', save = 'dotplot-mature-exp.pdf')


# In[24]:


scv.pl.velocity_embedding(f, basis = 'umap', color = 'louvain', save = 'velocity-embedding.pdf')


# In[25]:


scv.pl.velocity_embedding_stream(f, basis = 'umap', color = "louvain",save = 'velocity-embedding-stream-steady.pdf')


# In[26]:


mask = (f.obs["louvain"]  == "3") | (f.obs["louvain"]  == "1")
final = f[mask].copy()


# In[27]:


scv.pp.neighbors(final)
scv.pp.moments(final, n_pcs=30, n_neighbors=30)
scv.tl.velocity(final)
scv.tl.velocity_graph(final)
scv.tl.umap(final, n_components=2)
scv.tl.louvain(final)


# In[28]:


f


# In[29]:


scv.pl.velocity_embedding_stream(final, basis = 'umap', color = "louvain",save = 'velocity-embedding-stream-steady.pdf')


# In[30]:


cell = final.obs_names.tolist()
fdata = ddata[cell,:]


# In[31]:


c = final.obs['louvain']
fdata.obs['louvains'] = c


# In[32]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Gap43"]
scp.pl.violin(fdata, keys=immature_neuron_markers, groupby='louvains', save = 'immature_neurons2_exp.pdf')


# In[33]:


scp.pl.dotplot(fdata, immature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvains', save = 'dotplot-immature2-exp.pdf')


# In[34]:


#Expression of Immature neuron marker to assign clusters their stage
immature_neuron_markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43"]
scp.pl.violin(fdata, keys=immature_neuron_markers, groupby='louvains', save = 'immature_neurons3_exp.pdf')


# In[35]:


scp.pl.dotplot(fdata, immature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvains', save = 'dotplot-immature3-exp.pdf')


# In[36]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Omp"]
scp.pl.violin(fdata, keys=mature_neuron_markers, groupby='louvains', save = 'mature2-exp.pdf')


# In[37]:


scp.pl.dotplot(fdata, mature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvains', save = 'dotplot-mature2-exp.pdf')


# In[38]:


#Mature neuron markers to name the clusters on basis of stage
mature_neuron_markers = ["Myo6", "Stoml3",
                           "Gng13","Ech1","Dcx","Gnal", "Omp"]
scp.pl.violin(fdata, keys=mature_neuron_markers, groupby='louvains', save = 'mature3-exp.pdf')


# In[39]:


scp.pl.dotplot(fdata, mature_neuron_markers,dendrogram=True, cmap='Blues', groupby='louvains', save = 'dotplot-mature3-exp.pdf')


# #Final cluster names
# 

# In[40]:


final.obs['celltype'] = 'Immature cells'
final.obs['celltype'][final.obs['louvain']=='1'] = 'Mature cells'
final.obs['celltype'][final.obs['louvain']=='6'] = 'Transition cells'
final.obs['celltype'][final.obs['louvain']=='2'] = 'Mature cells'
final.obs['celltype'][final.obs['louvain']=='4'] = 'Mature cells'
final.obs['celltype'][final.obs['louvain']=='5'] = 'Mature cells'


# In[41]:


scv.pl.velocity_embedding_stream(final, basis = 'umap', color = 'celltype',save= 'velocity_embedding_stream.pdf')


# In[42]:


k = final.obs['celltype']
ddata.obs['celltype'] = k


# In[43]:


markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43", "Myo6", "Stoml3",
                           "Gng13","Ech1","Dcx","Gnal"]


# In[44]:


scp.pl.dotplot(ddata, markers,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot_all_markers2.pdf')


# EXTRACTING BARCODES

# In[45]:


f = final


# In[46]:


f


# In[47]:


trans = f[f.obs['celltype'] == 'Transition cells']
mature = f[f.obs['celltype'] == 'Mature cells']
immature = f[f.obs['celltype'] == 'Immature cells']


# In[48]:


trans_bar = pd.DataFrame(trans.obs.index.tolist())
mat_bar = pd.DataFrame(mature.obs.index.tolist())
imm_bar = pd.DataFrame(immature.obs.index.tolist())


# In[51]:


imm_bar.shape


# In[52]:


trans_bar.to_csv('transition_barcodes.csv')
mat_bar.to_csv('mature_barcodes.csv')
imm_bar.to_csv('immature_barcodeds')


# # CLUSTER SEGREGATION TO FIND THEIR IDENTITY IN R
# 
# 
# 

# IMMATURE

# In[53]:


clus = (f.obs["celltype"]  == "Immature cells")
clus1 = f[clus].copy()


# In[54]:


clus1_ind = clus1.obs_names


# In[55]:


clus1_list = clus1_ind.to_list()
clust1 = adata[clus1_list, :]


# Umap coordinates

# In[56]:


umap1 = clus1.obsm['X_umap']
umap1 = pd.DataFrame(umap1)
umap1.index = clus1.obs_names


# In[57]:


umap1.to_csv('umapimmat.csv')


# spliced of mature --- 1

# In[58]:


spl1 = clust1.layers['spliced']
spl1 = spl1.todense()
spl1 = pd.DataFrame(spl1)


# In[59]:


spl1.columns = clust1.var_names
spl1.index = clust1.obs_names


# In[60]:


spl1.to_csv('splimmat.csv')


# Unspliced of mature

# In[61]:


uspl1 = clust1.layers['unspliced']
uspl1 = uspl1.todense()
uspl1 = pd.DataFrame(uspl1)


# In[62]:


uspl1.columns = clust1.var_names
uspl1.index = clust1.obs_names


# In[63]:


uspl1.to_csv('usplimma.csv')


# TRANSITION

# In[64]:


#cluster0 TRANSITION
clus = (f.obs["celltype"]  == "Transition cells")
clus0 = f[clus].copy()


# In[65]:


clus0_ind = clus0.obs_names
clus0_list = clus0_ind.to_list()
#clus0_list = [sub[ : -2] for sub in clus0_ind]
clust0 = adata[clus0_list, :]


# umap

# In[66]:


umap0 = clus0.obsm['X_umap']
umap0 = pd.DataFrame(umap0)
umap0.index = clus0.obs_names


# In[67]:


umap0.to_csv('umaptrans.csv')


# spliced

# In[68]:


spl0 = clust0.layers['spliced']
spl0 = spl0.todense()
spl0 = pd.DataFrame(spl0)
spl0.columns = clust0.var_names
spl0.index = clust0.obs_names


# In[69]:


spl0.to_csv('spltrans.csv')


# In[70]:


uspl0 = clust0.layers['unspliced']
uspl0 = uspl0.todense()
uspl0 = pd.DataFrame(uspl0)


# In[71]:


uspl0.columns = clust0.var_names
uspl0.index = clust0.obs_names


# In[72]:


uspl0.to_csv('uspltrans.csv')


# #cluster mature

# In[73]:


clus = (f.obs["celltype"]  == "Mature cells")
clus4 = f[clus].copy()


# In[74]:


clus4_ind = clus4.obs_names
clus4_list = clus4_ind.to_list()
#clus4_list = [sub[ : -2] for sub in clus4_ind]
clust4 = adata[clus4_list, :]


# umap

# In[75]:


umap4 = clus4.obsm['X_umap']
umap4 = pd.DataFrame(umap4)
umap4.index = clus4.obs_names


# In[76]:


umap4.to_csv('umapmature.csv')


# spliced

# In[77]:


spl4 = clust4.layers['spliced']
spl4 = spl4.todense()
spl4 = pd.DataFrame(spl4)


# In[78]:


spl4.columns = clust4.var_names
spl4.index = clust4.obs_names


# In[79]:


spl4.to_csv('splmature.csv')


# unspliced

# In[80]:


uspl4 = clust4.layers['unspliced']
uspl4 = uspl4.todense()
uspl4 = pd.DataFrame(uspl4)


# In[81]:


uspl4.columns = clust4.var_names
uspl4.index = clust4.obs_names


# In[82]:


uspl4.to_csv('usplmature.csv')


# # Downstream analysis
# IMPORTANT GENES

# In[53]:


#file with all genes
ddata = f
ddata.write_loom("downstream.loom", write_obsm_varm=True)


# In[54]:


scv.tl.rank_velocity_genes(f, groupby='celltype', min_corr=.3)
df = scv.DataFrame(f.uns['rank_velocity_genes']['names'])
df.head()


# In[55]:


df.to_csv('rank_velocity_genes_wild.csv')


# VELOCITIES IN CYCLING PROGENITOR

# In[60]:


gene = f.var_names.tolist()


# In[81]:


scv.logging.print_versions()


# In[82]:


scv.utils.get_phase_marker_genes(f)


# In[94]:


scv.tl.score_genes_cell_cycle(f, s_genes = ['Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
     'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Mlf1ip', 'Hells', 'Rfc2',
     'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
     'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
     'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8'])
scv.pl.scatter(f, color_gradients=['S_score', 'G2M_score' ], perc=[5, 95], save = 'cell-score.pdf')


# SPEED AND COHERENCE

# In[62]:


f


# In[63]:


scv.tl.velocity_confidence(f)


# In[64]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(f, c=keys, cmap='coolwarm', perc=[5, 95], color = 'celltype', save = 'speed_and_coherence.pdf')


# In[65]:


df = f.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)


# In[66]:


df.to_csv('velocity_data.csv')


# VELOCITY GRAPH AND PSEUDOTIME

# In[67]:


scv.pl.velocity_graph(f, threshold=.1, color='celltype', save = 'Velocity_graph.pdf')


# In[68]:


x, y = scv.utils.get_cell_transitions(f, basis='umap', starting_cell= 100)
ax = scv.pl.velocity_graph(f, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(f, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save='descendants-ancestor-cell.pdf')


# In[69]:


scv.tl.velocity_pseudotime(f)
scv.pl.scatter(f, color='velocity_pseudotime', cmap='gnuplot', save = 'pseudotime.pdf')


# PAGA VELOCITY GRAPH

# In[70]:


f.uns['neighbors']['distances'] = f.obsp['distances']
f.uns['neighbors']['connectivities'] = f.obsp['connectivities']
scv.tl.paga(f, groups='celltype')
df = scv.get_df(f, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[71]:


scv.pl.paga(f, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save = 'PAGA.pdf')


# # Dynamical modeling

# In[72]:


scv.tl.recover_dynamics(f)


# In[73]:


scv.tl.velocity(f, mode='dynamical')
scv.tl.velocity_graph(f)


# In[74]:


scv.pl.velocity_embedding_stream(f, basis='umap', color = 'celltype', save = 'differential_velocity_stream.pdf')


# In[77]:


scv.tl.latent_time(f)
scv.pl.scatter(f, color='latent_time', color_map='gnuplot', size=80, save = 'Latent_time.pdf')


# CLUSTER SPECIFIC TOP LIKELIHOOD GENES

# In[80]:


scv.tl.rank_dynamical_genes(f, groupby='celltype')
df = scv.get_df(f, 'rank_dynamical_genes/names')
df.head(5)


# In[81]:


df.to_csv('rank_dynamical_genes.csv')


# # CELL RANK

# In[83]:


cr.tl.terminal_states(f, weight_connectivities=0.2,cluster_key='celltype')


# In[84]:


cr.pl.terminal_states(f, cluster_key='celltype', save = 'cellrank_Terminalstates.pdf')


# In[85]:


cr.tl.initial_states(f, cluster_key='celltype')
cr.pl.initial_states(f, discrete=True, save = 'cellrank_initialstates.pdf')


# In[86]:


cr.tl.lineages(f)
cr.pl.lineages(f, same_plot=False,rescale_color = [0,1], save = 'cellrank_lineages.pdf')


# In[87]:


scv.tl.latent_time(f, root_key= 'initial_state_probs', end_key='terminal_states_probs')


# In[88]:


scv.tl.paga(f, groups='celltype', root_key='initial_state_probs', end_key='terminal_states_probs',
            use_time_prior='velocity_pseudotime')


# In[90]:


cr.tl.lineage_drivers(f)


# In[91]:


root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]


# In[92]:


f.uns['iroot'] = root_idx


# In[93]:


scp.tl.dpt(f)


# In[94]:


# compue dpt, starting from CellRank defined root cell
scv.pl.scatter(f, color=['celltype', root_idx, 'latent_time', 'dpt_pseudotime'], fontsize=16,
               cmap='viridis', perc=[2, 98], colorbar=True, rescale_color=[0, 1],
               title=['louvain', 'root cell', 'latent time', 'dpt pseudotime'], save = 'cellrank_combinedgraph.pdf')


# In[95]:


model = cr.ul.models.GAM(f)


# ADVANCE

# In[97]:


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(f)


# In[98]:


vk.compute_transition_matrix()


# In[99]:


from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(f).compute_transition_matrix()


# In[100]:


combined_kernel = 0.8 * vk + 0.2 * ck
root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]
f.uns['iroot'] = root_idx


# In[101]:


from cellrank.tl.kernels import PalantirKernel
pk = PalantirKernel(f, time_key='dpt_pseudotime').compute_transition_matrix()
print(pk)


# In[102]:


from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)
print(g)


# In[105]:


g.compute_macrostates(n_states=2, cluster_key='celltype')
g.plot_macrostates(save = 'cellrank-macrostates.pdf')


# In[106]:


g.plot_macrostates(same_plot=False, rescale_color = [0,1], save = 'cellrank-macrostates-diff-plots.pdf')


# In[107]:


g.plot_macrostates(discrete=True, save = 'cellrank_macrostates-forward.pdf')


# In[108]:


g.set_terminal_states_from_macrostates()


# In[109]:


g.compute_absorption_probabilities()


# In[110]:


cr.pl.cluster_fates(f, cluster_key='celltype', mode='heatmap', save = 'cellrank-heatmap-clusterfate.pdf')


# # Differential kinetics

# In[111]:


scv.tl.differential_kinetic_test(f, groupby='celltype')


# In[115]:


scv.tl.recover_dynamics(f)


# In[116]:


top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(f, groupby='celltype')


# RECOMPUTING VELOCITIES

# In[119]:


scv.tl.velocity(f, diff_kinetics=True, groupby='celltype')
scv.tl.velocity_graph(f)


# In[120]:


scv.pl.velocity_embedding(f, dpi=120, arrow_size=2, arrow_length=2, save = 'differential-velocity.pdf')


# DGE

# In[121]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=10) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'differentialexp.pdf')


# In[123]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=5) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'diffexp5genes.pdf')


# In[125]:


im = ['Stmn1', 'Sox11', 'Gap43', 'Gng8', 'Vim']
scp.pl.dotplot(f, im,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-immature-diffexp.pdf')


# In[127]:


tr = ['Stmn4', 'Gng13', 'S100a5', 'Mfge8', 'Kirrel2']
scp.pl.dotplot(f, tr,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-trans-diffexp.pdf')


# In[128]:


m = ['Stoml3', 'Umodl1', 'Stom', 'Gnal', 'Chga']
scp.pl.dotplot(f, m,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-mature-diffexp.pdf')


# velocity length and confidence

# In[129]:


vl = f.obs['velocity_length']
vl.to_csv('velocity_length_wild.csv')


# In[130]:


vc = f.obs['velocity_confidence']
vc.to_csv('velocity_confi_wild.csv')

