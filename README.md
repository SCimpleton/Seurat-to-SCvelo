# Seurat-to-SCvelo
SCimple way of reproducing your Seurat UMAP with RNA velocities overlayed

```markdown
##NB make sure special characters are removed from cluster IDs##

#view the metadata to get a handle on the cell IDs seurat has created#
limb1@meta.data$cell<- rownames(limb1@meta.data)
meta<- limb1@meta.data
```
At this stage you will need to replace the numbers and underscores to make each cell the same as it appears in the velocyto loom output.
Velocyto gets round the duplicate cell names by adding a batch key from the 10x channel. For some reason it then seems to add an appended 'x'
so swap appended characters for the corresponding batch (hopefully this was metadata in your object)

```markdown
limb1@meta.data$cell<- gsub('1_1_1_1', '1-lane1', limb1@meta.data$cell)
limb1@meta.data$cell<- gsub('1_2_1_1', '1-lane2', limb1@meta.data$cell)
limb1@meta.data$cell<- gsub('1_1_1', '1-lane3', limb1@meta.data$cell)
limb1@meta.data$cell<- gsub('1_2_1', '1-lane4', limb1@meta.data$cell)
limb1@meta.data$cell<- gsub('1_1', '1-lane5', limb1@meta.data$cell)
limb1@meta.data$cell<- gsub('1_2', '1-lane6', limb1@meta.data$cell)

rownames(limb1@meta.data)<- limb1@meta.data$cell

#then tweak the formatting to match velocyto
sample_obs<- limb1@meta.data$cell
sample_obs<- as.data.frame(sample_obs)
sample_obs<- separate(sample_obs, sample_obs, sep = "-1-",into = c("barcode", "study"))
sample_obs<- unite(sample_obs, study, barcode, col = "x", sep = ":")
sample_obs$x <- paste0(sample_obs$x, "x" , sep="")

#write
write.csv(sample_obs, file = "cellID_obs_seurat.csv", row.names = FALSE)

#Now repeat this process for UMAP embeddings
umap<- (Embeddings(limb1, reduction = "umap"))
rownames<- rownames(limb1@meta.data)
umap<- as.data.frame(umap, rownames=F)
umap$cell<- rownames
umap<- separate(umap,cell , sep = "-1-",into = c("barcode", "study"))
umap<- unite(umap, study, barcode, col = "cell", sep = ":")
umap$cell <- paste0(umap$cell, "x" , sep="")
write.csv(umap, file = "cell_embeddings_seurat.csv", row.names = F)

#and the same for cluster info
meta<- limb1@meta.data

#choose the annotation column here
clusters<- meta[,5]
clusters<- as.data.frame(clusters)
clusters$cell<- rownames
clusters<- separate(clusters,cell , sep = "-1-",into = c("barcode", "study"))
clusters<- unite(clusters, study, barcode, col = "cell", sep = ":")
clusters$cell <- paste0(clusters$cell, "x" , sep="")

#write this initial table as the clusters
write.csv(clusters, file = "cluster_annot_seurat.csv", row.names = F)

#finally, convert the cluster names into the colours you have used for each in the seurat umap plots
#now convert cluster to colour
cols= pals::cols25(n=25)

clusters$clusters <- gsub("celltypeA", "#1F78C8", clusters$clusters)
clusters$clusters <- gsub("celltypeB", "#ff0000", clusters$clusters)
clusters$clusters <- gsub("celltypeC", "#33a02c", clusters$clusters)
#etc etc

#write the file
write.csv(clusters, file = "clusters_obs_seurat.csv", row.names = F)
```


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

#now subset by removing irrelevant cells form the loom file from velocyto
sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["x"])]

#check all removed appropriately ie remaining cells mateches those in the seurat object
sample_one.obs.index

#convert index to dataframe
sample_one_index = pd.DataFrame(sample_one.obs.index)

#rename column to cellID
sample_one_index = sample_one_index.rename(columns = {0:'Cell ID'})

#do the same for umap data, then order appropriately
umap = umap.rename(columns = {'cell':'Cell ID'})
umap_ordered = sample_one_index.merge(umap, on = "Cell ID")

#inspect
umap_ordered

#remove the cellID column - check the location is 1 above - then add as an observation in the anndata
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values

#repeat for clusters, but as an unstructured annotation
cell_clusters = cell_clusters.rename(columns = {'cell':'Cell ID'})
cell_clusters_ordered = sample_one_index.merge(cell_clusters, on = "Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
sample_one.uns['Clusters'] = cell_clusters_ordered.values

#compute
scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)

#plot - lots of options here but basic one can be made like this
scv.pl.velocity_embedding_grid(sample_one, basis = 'X_umap', color = sample_one.uns['Clusters'],arrow_size=2,arrow_color="black", scale=2,figsize=(6.5,7), arrow_length=5)

#save if happy
scv.pl.velocity_embedding_stream(sample_one, basis = 'X_umap', color = sample_one.uns['Clusters'], min_mass=0, figsize=(6.5,7), save=".pdf")

#other functions can be carried out from here such as inspecting gene dynamics- see scvelo documentation
