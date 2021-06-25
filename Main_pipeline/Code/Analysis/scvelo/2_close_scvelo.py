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
os.chdir('/home/sidrah19220/nmd/mouse/analysis/scvelo/close')


# In[4]:


scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
scv.set_figure_params()


# In[5]:


file = "close.loom"
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


mask = (f.obs["louvain"]  == "8") | (f.obs["louvain"]  == "2")
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


final.obs['celltype'] = 'Mature cells'
final.obs['celltype'][final.obs['louvain']=='3'] = 'Immature cells'
final.obs['celltype'][final.obs['louvain']=='4'] = 'Transition cells'


# In[41]:


scv.pl.velocity_embedding_stream(final, basis = 'umap', color = 'celltype',save= 'velocity_embedding_stream.pdf')


# In[42]:


k = final.obs['celltype']
ddata.obs['celltype'] = k


# In[43]:


markers = ["Ncam1","Lhx2","Gng8", "Rcor1","Ctxn3","Gap43", "Myo6", "Stoml3",
                           "Gng13","Ech1","Dcx","Gnal", "Omp"]


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


# In[49]:


trans_bar.shape


# In[50]:


mat_bar.shape


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

# In[69]:


spl0 = clust0.layers['spliced']
spl0 = spl0.todense()
spl0 = pd.DataFrame(spl0)
spl0.columns = clust0.var_names
spl0.index = clust0.obs_names


# In[70]:


spl0.to_csv('spltrans.csv')


# In[71]:


uspl0 = clust0.layers['unspliced']
uspl0 = uspl0.todense()
uspl0 = pd.DataFrame(uspl0)


# In[72]:


uspl0.columns = clust0.var_names
uspl0.index = clust0.obs_names


# In[73]:


uspl0.to_csv('uspltrans.csv')


# #cluster mature

# In[74]:


clus = (f.obs["celltype"]  == "Mature cells")
clus4 = f[clus].copy()


# In[75]:


clus4_ind = clus4.obs_names
clus4_list = clus4_ind.to_list()
#clus4_list = [sub[ : -2] for sub in clus4_ind]
clust4 = adata[clus4_list, :]


# umap

# In[76]:


umap4 = clus4.obsm['X_umap']
umap4 = pd.DataFrame(umap4)
umap4.index = clus4.obs_names


# In[77]:


umap4.to_csv('umapmature.csv')


# spliced

# In[78]:


spl4 = clust4.layers['spliced']
spl4 = spl4.todense()
spl4 = pd.DataFrame(spl4)


# In[79]:


spl4.columns = clust4.var_names
spl4.index = clust4.obs_names


# In[80]:


spl4.to_csv('splmature.csv')


# unspliced

# In[81]:


uspl4 = clust4.layers['unspliced']
uspl4 = uspl4.todense()
uspl4 = pd.DataFrame(uspl4)


# In[82]:


uspl4.columns = clust4.var_names
uspl4.index = clust4.obs_names


# In[83]:


uspl4.to_csv('usplmature.csv')


# # Downstream analysis
# IMPORTANT GENES

# In[52]:


#file with all genes
ddata = f
ddata.write_loom("downstream.loom", write_obsm_varm=True)


# In[53]:


scv.tl.rank_velocity_genes(f, groupby='celltype', min_corr=.3)
df = scv.DataFrame(f.uns['rank_velocity_genes']['names'])
df.head()


# In[54]:


df.to_csv('rank_velocity_genes_wild.csv')


# In[55]:


#Expressions of velocity genes
kwargs = dict(frameon=False, size=10, linewidth=1.5)

scv.pl.scatter(f, df['Immature cells'][:5], ylabel='Immature cells', **kwargs, save = 'scatterplot_Immature_cells.pdf')
scv.pl.scatter(f, df['Mature cells'][:5], ylabel='Mature cells', **kwargs, save = 'scatterplot_mature_cells.pdf')
scv.pl.scatter(f, df['Mature cells'][:5], ylabel='Transition cells', **kwargs, save = 'scatterplot_transition_cells.pdf')


# VELOCITIES IN CYCLING PROGENITOR

# In[62]:


scv.tl.score_genes_cell_cycle(f)
scv.pl.scatter(f, color_gradients=['S_score', 'G2M_score' ], perc=[5, 95], save = 'cell-score.pdf')


# In[90]:


scv.pl.velocity(f, ['Gng8', 'Gap43', 'Ctxn3'], add_outline=True, save= 'Velocitygraph_Immature_cell.pdf')
scv.pl.velocity(f, ['Gnal', 'Dcx', 'Camp'], add_outline = True, save= 'Velocitygraph_mature_cell.pdf')


# SPEED AND COHERENCE

# In[91]:


f


# In[92]:


scv.tl.velocity_confidence(f)


# In[93]:


keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(f, c=keys, cmap='coolwarm', perc=[5, 95], color = 'celltype', save = 'speed_and_coherence.pdf')


# In[94]:


df = f.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)


# In[95]:


df.to_csv('velocity_data.csv')


# VELOCITY GRAPH AND PSEUDOTIME

# In[96]:


scv.pl.velocity_graph(f, threshold=.1, color='celltype', save = 'Velocity_graph.pdf')


# In[97]:


x, y = scv.utils.get_cell_transitions(f, basis='umap', starting_cell= 100)
ax = scv.pl.velocity_graph(f, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(f, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save='descendants-ancestor-cell.pdf')


# In[98]:


scv.tl.velocity_pseudotime(f)
scv.pl.scatter(f, color='velocity_pseudotime', cmap='gnuplot', save = 'pseudotime.pdf')


# PAGA VELOCITY GRAPH

# In[99]:


f.uns['neighbors']['distances'] = f.obsp['distances']
f.uns['neighbors']['connectivities'] = f.obsp['connectivities']
scv.tl.paga(f, groups='celltype')
df = scv.get_df(f, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[100]:


scv.pl.paga(f, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save = 'PAGA.pdf')


# # Dynamical modeling

# In[101]:


scv.tl.recover_dynamics(f)


# In[102]:


scv.tl.velocity(f, mode='dynamical')
scv.tl.velocity_graph(f)


# In[103]:


scv.pl.velocity_embedding_stream(f, basis='umap', color = 'celltype', save = 'differential_velocity_stream.pdf')


# In[104]:


df = f.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate',xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)
    plt.savefig('trans-spl-deg.pdf')

scv.get_df(f, 'fit*', dropna=True).head()


# In[105]:


df.to_csv('trans-spli-degr-wild.csv')


# In[106]:


scv.tl.latent_time(f)
scv.pl.scatter(f, color='latent_time', color_map='gnuplot', size=80, save = 'Latent_time.pdf')


# In[107]:


top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.pl.heatmap(f, var_names=top_genes, sortby='latent_time', n_convolve=100, save = 'Toplikelihoodgenes_heatmap.pdf')


# TOP LIKELIHOOD GENES - DRIVER GENES

# In[108]:


top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(f, basis=top_genes[:15], ncols=5, frameon=False, color = 'celltype',save = 'toplikelihood-15genes_scatterplot')


# CLUSTER SPECIFIC TOP LIKELIHOOD GENES

# In[109]:


scv.tl.rank_dynamical_genes(f, groupby='celltype')
df = scv.get_df(f, 'rank_dynamical_genes/names')
df.head(5)


# In[110]:


df.to_csv('rank_dynamical_genes.csv')


# In[111]:


for cluster in ['Immature cells', 'Mature cells']:
    scv.pl.scatter(f, df[cluster][:5], ylabel=cluster, frameon=False, color = 'celltype', save= 'cluster-specific_likelihood_genes{x}'+cluster)


# # CELL RANK

# In[112]:


cr.tl.terminal_states(f, weight_connectivities=0.2,cluster_key='celltype')


# In[113]:


cr.pl.terminal_states(f, cluster_key='celltype', save = 'cellrank_Terminalstates.pdf')


# In[114]:


cr.tl.initial_states(f, cluster_key='celltype')
cr.pl.initial_states(f, discrete=True, save = 'cellrank_initialstates.pdf')


# In[115]:


cr.tl.lineages(f)
cr.pl.lineages(f, same_plot=False,rescale_color = [0,1], save = 'cellrank_lineages.pdf')


# In[116]:


scv.tl.latent_time(f, root_key= 'initial_state_probs', end_key='terminal_states_probs')


# In[117]:


scv.tl.paga(f, groups='celltype', root_key='initial_state_probs', end_key='terminal_states_probs',
            use_time_prior='velocity_pseudotime')


# In[118]:


cr.pl.cluster_fates(f, mode="paga_pie", cluster_key="celltype", basis='umap', 
                    legend_kwargs={'loc': 'top right out'}, legend_loc='top left out', 
                    node_size_scale=5, edge_width_scale=1, max_edge_width=4, title='directed PAGA', save = 'Cellrank-Differential_PAGA.pdf')


# In[119]:


cr.tl.lineage_drivers(f)


# In[120]:


root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]


# In[121]:


f.uns['iroot'] = root_idx


# In[122]:


scp.tl.dpt(f)


# In[123]:


# compue dpt, starting from CellRank defined root cell
scv.pl.scatter(f, color=['celltype', root_idx, 'latent_time', 'dpt_pseudotime'], fontsize=16,
               cmap='viridis', perc=[2, 98], colorbar=True, rescale_color=[0, 1],
               title=['louvain', 'root cell', 'latent time', 'dpt pseudotime'], save = 'cellrank_combinedgraph.pdf')


# In[124]:


model = cr.ul.models.GAM(f)


# In[125]:


cr.pl.gene_trends(f, model=model, data_key='X',
                  genes=['Gng8', 'Gap43'], ncols=2,
                  time_key='latent_time', same_plot=True, hide_cells=True,
                  figsize=(15, 4), n_test_points=200)


# ADVANCE

# In[126]:


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(f)


# In[127]:


vk.compute_transition_matrix()


# In[128]:


from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(f).compute_transition_matrix()


# In[129]:


combined_kernel = 0.8 * vk + 0.2 * ck
root_idx = np.where(f.obs['celltype'] == 'Immature cells')[0][0]
f.uns['iroot'] = root_idx


# In[130]:


from cellrank.tl.kernels import PalantirKernel
pk = PalantirKernel(f, time_key='dpt_pseudotime').compute_transition_matrix()
print(pk)


# In[131]:


from cellrank.tl.estimators import GPCCA
g = GPCCA(combined_kernel)
print(g)


# In[132]:


g.compute_schur(n_components=20)
g.plot_spectrum(save='cellrank-eigenvalues.pdf')


# In[133]:


g.plot_schur(use=1, rescale_color = [0,1], save='cellrank-schurvector.pdf')


# In[134]:


g.compute_macrostates(n_states=2, cluster_key='celltype')
g.plot_macrostates(save = 'cellrank-macrostates.pdf')


# In[135]:


g.plot_macrostates(same_plot=False, rescale_color = [0,1], save = 'cellrank-macrostates-diff-plots.pdf')


# In[136]:


g.plot_macrostates(discrete=True, save = 'cellrank_macrostates-forward.pdf')


# In[137]:


g.set_terminal_states_from_macrostates()


# In[138]:


g.compute_absorption_probabilities()


# In[139]:


cr.pl.cluster_fates(f, cluster_key='celltype', mode='heatmap', save = 'cellrank-heatmap-clusterfate.pdf')


# # Differential kinetics

# In[140]:


scv.tl.differential_kinetic_test(f, groupby='celltype')


# In[141]:


var_names = ['Gap43', 'Gng8']


# In[142]:


kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(f, basis=var_names, add_outline='fit_diff_kinetics', **kwargs)


# In[143]:


diff_clusters=list(f[:, var_names].var['fit_diff_kinetics'])
scv.pl.scatter(f, legend_loc='right', size=60, title='diff kinetics',
               add_outline=diff_clusters)


# TOP LIKELIHOOD genes

# In[144]:


scv.tl.recover_dynamics(f)


# In[145]:


top_genes = f.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(f, groupby='celltype')


# In[146]:


scv.pl.scatter(f, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics', save = 'differential-kinetics-15genes.pdf')


# In[147]:


scv.pl.scatter(f, basis=top_genes[15:30], ncols=5, add_outline='fit_diff_kinetics', save = 'differential-kinetics-15-30genes.pdf')


# RECOMPUTING VELOCITIES

# In[148]:


scv.tl.velocity(f, diff_kinetics=True, groupby='celltype')
scv.tl.velocity_graph(f)


# In[149]:


scv.pl.velocity_embedding(f, dpi=120, arrow_size=2, arrow_length=2, save = 'differential-velocity.pdf')


# DGE

# In[150]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=10) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'differentialexp.pdf')


# In[151]:


scp.tl.rank_genes_groups(f, groupby='celltype', use_raw=True, 
                        method='t-test_overestim_var', n_genes=5) # compute differential expression
scp.pl.rank_genes_groups_tracksplot(f, groupby='celltype', save = 'diffexp5genes.pdf')


# In[153]:


im = ['Stmn1', 'Sox11', 'Gap43', 'Gng8', 'Vim']
scp.pl.dotplot(f, im,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-immature-diffexp.pdf')


# In[154]:


tr = ['Insm1', 'Gng8', 'Stmn2', 'Stmn1', 'Gap43']
scp.pl.dotplot(f, tr,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-trans-diffexp.pdf')


# In[155]:


m = ['Stoml3', 'Umodl1', 'Dcdc2a', 'Nsg1', 'Gng13']
scp.pl.dotplot(f, m,dendrogram=True, cmap='Blues', groupby='celltype', save = 'dotplot-mature-diffexp.pdf')


# velocity length and confidence

# In[156]:


vl = f.obs['velocity_length']
vl.to_csv('velocity_length_wild.csv')


# In[157]:


vc = f.obs['velocity_confidence']
vc.to_csv('velocity_confi_wild.csv')

