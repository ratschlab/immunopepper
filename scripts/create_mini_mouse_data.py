import pickle
import numpy as np
import h5py
## create mini-version splicegraph
a = pickle.load(open('spladder_out/spladder/genes_graph_conf3.merge_graphs.pickle','rb'),encoding='latin1')
choose_gene_id = [19,33]
strain_id = [0,1] # ENCSR000BZG ERR2130621
new_graph_tuple = (a[0][choose_gene_id],a[1])
f = open('ImmunoPepper_usecase.pickle','wb')
pickle.dump(new_graph_tuple,f)
f.close()

## create mini-version count file
h5 = h5py.File('spladder_out/spladder/genes_graph_conf3.merge_graphs.count.hdf5','r')
newcount = h5py.File('ImmunoPepper_usecase.count.hdf5','w')
gene_ids_edges = [i for i,item in enumerate(h5['gene_ids_edges'][:1000]) if item in choose_gene_id ]
gene_ids_segs = [i for i,item in enumerate(h5['gene_ids_segs'][:1000]) if item in choose_gene_id ]

break_id = np.where(np.diff(gene_ids_edges)>1)[0][0]
new_ids_edges = np.reshape(gene_ids_edges,(len(gene_ids_edges),1))
new_ids_edges[:break_id+1] = 0
new_ids_edges[break_id+1:] = 1

break_id = np.where(np.diff(gene_ids_segs)>1)[0][0]
new_ids_segs = np.reshape(gene_ids_segs,(len(gene_ids_segs),1))
new_ids_segs[:break_id+1] = 0
new_ids_segs[break_id+1:] = 1

ind = gene_ids_edges
newcount.create_dataset('edge_idx',data=h5['edge_idx'][ind])
newcount.create_dataset('edges',data=h5['edges'][ind][:,strain_id])
newcount.create_dataset('gene_ids_edges',data=new_ids_edges)
newcount.create_dataset('gene_names',data=h5['gene_names'][choose_gene_id])

ind = gene_ids_segs
newcount.create_dataset('gene_ids_segs',data=new_ids_segs)
newcount.create_dataset('segments',data=h5['segments'][ind][:,strain_id])
newcount.create_dataset('seg_pos',data=h5['seg_pos'][ind][:,strain_id])
newcount.create_dataset('strains',data=h5['strains'][strain_id])
newcount.close()
