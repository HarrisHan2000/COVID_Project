#!/usr/bin/env python
# coding: utf-8

# In[41]:


# Import essential libraries
import os as os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white')


# # Read datasets, QC, and integrate them

# In[42]:


# Unify Vars
varnames = np.load(r'/home/aimmunelab/lab_members/Harris/COVID_Project/Datasets_selected/vars.npy',
                  allow_pickle=True)
def UnifyVars(adata, var = varnames):
    adata = adata[:,varnames]
    return adata


# In[43]:


# Set parental work directory
parental = '/home/aimmunelab/lab_members/Harris/COVID_Project/Datasets_selected'
os.chdir(parental)


# In[44]:


def BatchCategorization(groupname,samplelist, batchlist):
    if len(batchlist)==0:
        j = 0
    else:
        j = len(pd.Series(batchlist).astype('category').value_counts())
    for i in samplelist:
        batchlist+=[groupname+' '+str(j)]*len(i)
        j+=1


# In[45]:


# Define a basic QC metric function
def BasicQC(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    adata = adata[adata.obs.n_genes_by_counts > 750, :]
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[46]:


def Normalize(anndata):
    sc.pp.normalize_total(anndata, target_sum=1e4)
    sc.pp.log1p(anndata)
    anndata = anndata


# In[47]:


def Concatenate(anndatalist):
    concatenated = anndatalist[0].concatenate(anndatalist[1:len(anndatalist)])
    return concatenated


# In[48]:


# Prepare batch lists
healthy_batch = []
moderate_batch = []
severe_batch = []


# In[115]:


os.listdir(r'GSE163668_grouped/Severe')


# In[50]:


# Read GSE150728
os.chdir('GSE150728_grouped')
cata = os.listdir('Healthy')
healthy_group_GSE150728 = []
for i in cata:
    print('Reading '+i)
    sample = sc.read('Healthy/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    healthy_group_GSE150728.append(sample)
BatchCategorization('Healthy', healthy_group_GSE150728, healthy_batch)
healthy_group_GSE150728c = Concatenate(healthy_group_GSE150728)
healthy_group_GSE150728c.write('healthy_group_GSE150728c_raw.h5ad')
Normalize(healthy_group_GSE150728c)
healthy_group_GSE150728c.write('healthy_group_GSE150728c.h5ad')
# Severe Group
cata = os.listdir('Severe')
severe_group_GSE150728 = []
for i in cata:
    print('Reading'+i)
    sample = sc.read('Severe/'+i)
    BasicQC(sample)
    severe_group_GSE150728.append(sample)
BatchCategorization('Severe',severe_group_GSE150728,severe_batch)
severe_group_GSE150728c = Concatenate(severe_group_GSE150728)
severe_group_GSE150728c.write('severe_group_GSE150728c_raw.h5ad')
Normalize(severe_group_GSE150728c)
severe_group_GSE150728c.write('severe_group_GSE150728c.h5ad')
os.chdir(parental)


# In[51]:


# Read GSE166489
os.chdir('GSE166489_grouped')
# Healthy Group
cata = os.listdir('Healthy')
healthy_group_GSE166489 = []
for i in cata:
    sample = sc.read_10x_mtx('Healthy/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    healthy_group_GSE166489.append(sample)
BatchCategorization('Healthy', healthy_group_GSE166489, healthy_batch)
healthy_group_GSE166489c = Concatenate(healthy_group_GSE166489)
healthy_group_GSE166489c.write('healthy_group_GSE166489c_raw.h5ad')
Normalize(healthy_group_GSE166489c)
healthy_group_GSE166489c.write('healthy_group_GSE166489c.h5ad')
# Moderate Group
cata = os.listdir('Moderate')
moderate_group_GSE166489 = []
for i in cata:
    sample = sc.read_10x_mtx('Moderate/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    moderate_group_GSE166489.append(sample)
BatchCategorization('Moderate',moderate_group_GSE166489, moderate_batch)
moderate_group_GSE166489c = Concatenate(moderate_group_GSE166489)
moderate_group_GSE166489c.write('moderate_group_GSE166489c_raw.h5ad')
Normalize(moderate_group_GSE166489c)
moderate_group_GSE166489c.write('moderate_group_GSE166489c.h5ad')
# Severe Group
cata = os.listdir('Severe')
severe_group_GSE166489 = []
for i in cata:
    sample = sc.read_10x_mtx('Severe/'+i)
    sample = UnifyVars(sample)
    severe_group_GSE166489.append(sample)
    BasicQC(sample)
BatchCategorization('Severe',severe_group_GSE166489,severe_batch)
severe_group_GSE166489c = Concatenate(severe_group_GSE166489)
severe_group_GSE166489c.write('severe_group_GSE166489c_raw.h5ad')
Normalize(severe_group_GSE166489c)
severe_group_GSE166489c.write('severe_group_GSE166489c.h5ad')
os.chdir(parental)


# In[52]:


# Read GSE163668
os.chdir('GSE163668_grouped')
# Healthy Group
healthy_group_GSE163668 = []
for i in os.listdir('Healthy'):
    sample = sc.read_10x_mtx(r'Healthy/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    healthy_group_GSE163668.append(sample)
BatchCategorization('Healthy',healthy_group_GSE163668,healthy_batch)
healthy_group_GSE163668c = Concatenate(healthy_group_GSE163668)
Normalize(healthy_group_GSE163668c)
# Moderate Group
moderate_group_GSE163668 = []
cata = []
for i in os.listdir('Moderate'):
    if '_' not in i:
        cata.append(i)
for i in cata:
    sample = sc.read_10x_mtx(r'Moderate/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    moderate_group_GSE163668.append(sample)
BatchCategorization('Moderate', moderate_group_GSE163668,moderate_batch)
moderate_group_GSE163668c = Concatenate(moderate_group_GSE163668)
Normalize(moderate_group_GSE163668c)
# Severe Group
severe_group_GSE163668 = []
cata = []
for i in os.listdir('Severe'):
    if '_' not in i:
        cata.append(i)
for i in cata:
    sample = sc.read_10x_mtx(r'Severe/'+i)
    sample = UnifyVars(sample)
    BasicQC(sample)
    severe_group_GSE163668.append(sample)
BatchCategorization('Severe', severe_group_GSE163668, severe_batch)
severe_group_GSE163668c = Concatenate(severe_group_GSE163668)
Normalize(severe_group_GSE163668c)
os.chdir(parental)


# In[53]:


# Save batch lists
pd.DataFrame(healthy_batch).to_csv('healthy_batch.csv')
pd.DataFrame(moderate_batch).to_csv('moderate_batch.csv')
pd.DataFrame(severe_batch).to_csv('severe_batch.csv')


# In[54]:


healthy_batch = pd.read_csv('healthy_batch.csv')['0'].astype('category')
moderate_batch = pd.read_csv('moderate_batch.csv')['0'].astype('category')
severe_batch = pd.read_csv('severe_batch.csv')['0'].astype('category')


# In[55]:


# Combine batch lists and save
aggregate_batch = healthy_batch.append(moderate_batch)
aggregate_batch = aggregate_batch.append(severe_batch)
aggregate_batch = aggregate_batch.astype('category')
pd.DataFrame(aggregate_batch).to_csv('aggregate_batch.csv')


# In[56]:


# Transform lists into Series format with dytpe of category, which accords to original concatenated Anndatas
def Transform2CateSeries(batchlist):
    batchlist = pd.Series(batchlist).astype('category')
    return batchlist

aggregate_batch = Transform2CateSeries(aggregate_batch)
healthy_batch = Transform2CateSeries(healthy_batch)
moderate_batch = Transform2CateSeries(moderate_batch)
severe_batch = Transform2CateSeries(severe_batch)


# In[57]:


# Actually replace the batch item
def ReplaceBatch(batchlist,anndata):
    batchlist.index = anndata.obs['batch'].index
    anndata.obs['batch'] = batchlist


# In[58]:


healthy_group = healthy_group_GSE166489c.concatenate([healthy_group_GSE163668c, healthy_group_GSE150728c])
ReplaceBatch(healthy_batch,healthy_group)
moderate_group = moderate_group_GSE166489c.concatenate(moderate_group_GSE163668c)
ReplaceBatch(moderate_batch, moderate_group)
severe_group = severe_group_GSE166489c.concatenate([severe_group_GSE163668c, severe_group_GSE150728c])
ReplaceBatch(severe_batch,severe_group)


# In[59]:


healthy_group.write('healthy_group.h5ad')
moderate_group.write('moderate_group.h5ad')
severe_group.write('severe_group.h5ad')


# In[60]:


aggregate_batch = pd.Series()
for i in [healthy_batch, moderate_batch, severe_batch]:
    aggregate_batch = aggregate_batch.append(i)
aggregate_batch.astype('category')


# In[61]:


aggregate = healthy_group.concatenate([moderate_group,severe_group])
ReplaceBatch(aggregate_batch,aggregate)
aggregate.write('aggregate.h5ad')


# # Analyze the aggregate file

# In[67]:


aggregate_copy = sc.read('aggregate.h5ad')


# In[68]:


aggregate_copy.obs['batch']


# In[69]:


# Define a PCA-Harmony function
def PCAHarmony(adata):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='CST3')
    sc.external.pp.harmony_integrate(adata, 'batch')
# Define a neighborhood calculating and leiden clustering function
def Clustering(adata, pcs):
    sc.pp.neighbors(adata,n_pcs=pcs,use_rep='X_pca_harmony' )
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color='leiden')


# In[70]:


PCAHarmony(aggregate_copy)


# In[73]:


Clustering(aggregate_copy,10)


# In[74]:


sc.pl.umap(aggregate_copy,color='batch')


# In[75]:


aggregate_copy.write('aggregate_clustered.h5ad')


# In[110]:


len(aggregate_copy.obs['batch'].value_counts())


# # Quality control by sample volume (threshold = 2000)

# In[76]:


aggregate_copy = sc.read('aggregate_clustered.h5ad')
aggregate_copy


# In[80]:


sample_volume = []
sample_id = []
for i in aggregate_copy.obs['batch'].value_counts().index:
    print(i, len(aggregate_copy[aggregate_copy.obs['batch']==i]))
    sample_volume.append(len(aggregate_copy[aggregate_copy.obs['batch']==i]))
    sample_id.append(i)
sample_volume_df = pd.DataFrame({'Sample_id':sample_id,
              'Sample_volume':sample_volume})
sample_volume_df
sample_volume_df.to_csv('sample_volume_df.csv')


# In[81]:


# Keep samples contain more than 2000 cells
sample_volume_df[sample_volume_df['Sample_volume']>2000]


# In[96]:


aggregate_sample_qc = aggregate_copy
for i in sample_volume_df[sample_volume_df['Sample_volume']<2000]['Sample_id']:
    aggregate_sample_qc = aggregate_sample_qc[aggregate_sample_qc.obs['batch']!=str(i)]
aggregate_sample_qc


# In[83]:


sc.pl.umap(aggregate_sample_qc, color='leiden')


# # Quality control by cluster volume (threshold = 2000)

# In[85]:


aggregate_copy.obs['leiden'].value_counts()


# In[99]:


bad_cluster = np.arange(18,34)
for i in bad_cluster:
    aggregate_sample_qc = aggregate_sample_qc[aggregate_sample_qc.obs['leiden']!=str(i)]


# In[103]:


len(aggregate_sample_qc)


# In[100]:


sc.pl.umap(aggregate_sample_qc,color='leiden')


# In[91]:


aggregate_sample_qc.write('aggregate_clustered_QC_sample.h5ad')

