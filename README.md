# Seurat-to-SCvelo
SCimple way of reproducing your Seurat UMAP with RNA velocities overlayed






## Now switch to python
you will need to install the requisite packages if you havent already
```markdown
pip install -U scvelo
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

# Load each of the loom files generated from your single cell data
sample_one = anndata.read_loom("sample_1.loom")
sample_two = anndata.read_loom("sample_2.loom") 
sample_three = anndata.read_loom("sample_3.loom") 
sample_four = anndata.read_loom("sample_4.loom")
sample_five = anndata.read_loom("sample_5.loom")
sample_six = anndata.read_loom("sample_6.loom")
sample_seven = anndata.read_loom("sample_7.loom")
sample_eight = anndata.read_loom("sample_8.loom")
sample_nine = anndata.read_loom("sample_9.loom")
sample_ten = anndata.read_loom("sample_10.loom")

# merge them, telling the machine to ignore the non-unique indices. NB this final command varies depending on package version, so if it throws an error check that out
sample_one = sample_one.concatenate(sample_two, sample_three, sample_four, sample_five, sample_six, sample_seven, sample_eight, sample_nine, sample_ten,index_unique=None)

# now read in the files you made in R
sample_obs = pd.read_csv("cellID_obs_seurat.csv")
umap = pd.read_csv("cell_embeddings_seurat.csv")
cell_clusters = pd.read_csv("clusters_obs_seurat.csv")
cell_labels = pd.read_csv("cluster_annot_seurat.csv")
```

At this stage, it is worth a quick sense check that the loom file cells look the same as the cells you exported from seurat
so that the various bits of data can be correctly assigned to each cell, and the unwanted cells can be removed

```markdown
sample_one.obs.index
sample_obs
'''
