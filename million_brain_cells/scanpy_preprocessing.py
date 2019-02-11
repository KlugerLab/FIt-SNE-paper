import numpy
import matplotlib
#matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import scanpy.api as sc #Using scanpy 1.2.2

sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving

# Preprocess and find the Louvain clusters
filename = 'data/1M_neurons_filtered_gene_bc_matrices_h5.h5'
adata = sc.read_10x_h5(filename)
sc.pp.recipe_zheng17(adata) 
numpy.savetxt('data/processed_neurons.csv',adata.X, delimiter=',')


sc.pp.neighbors(adata,method='gauss')
sc.tl.louvain(adata)
f = open('data/louvain_gaussian.csv', 'w')
f.write("\n".join (adata.obs['louvain']))


# Write processed matrix
numpy.save('data/processed_neurons.npy',adata.X)


# We also need the raw  matrix so we can print out various genes
adata_ = sc.read_10x_h5(filename)
adata_.shape

gns2 = adata_.var_names.tolist()
gns2 = [x.lower() for x in gns2]


gois = ['sncg', 'slc17a8', 'spp1', 'col15a1'];
for goi in gois:
    idx = gns2.index(goi)
    numpy.savetxt('data/genes/%s.csv' % goi,adata_.X[:,idx].todense(), delimiter=',', fmt='%d')

