from base6 import *
for s in sections:
 au=js(N/s/'metadata/input_audit.json');sy=genedata(au)[1]
 gs=[direct,['FOXP3'],core[1:],epi,fib,mye]
 mapper=sparse.csr_matrix(np.array([np.isin(sy,gg) for gg in gs],np.int32).T);parts=[z@mapper for a,b,z in slabs(au)]
 sparse.save_npz(R/s/'observed_panel.npz',sparse.vstack(parts,format='csr'))
